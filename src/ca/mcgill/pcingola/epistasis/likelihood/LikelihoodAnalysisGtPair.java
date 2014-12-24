package ca.mcgill.pcingola.epistasis.likelihood;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.probablility.FisherExactTest;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.Genotype;
import ca.mcgill.pcingola.epistasis.gwas.GwasResult;
import ca.mcgill.pcingola.regression.LogisticRegression;
import ca.mcgill.pcingola.regression.LogisticRegressionIrwls;

/**
 * Logistic regression log-likelihood analysis of 2 VCF entries + phenotype data
 *
 * @author pcingola
 */
public class LikelihoodAnalysisGtPair extends LikelihoodAnalysisGt {

	public static final int MIN_SHARED_VARIANTS = 5;
	public static final double EPSILON = 1e-6;

	ArrayList<String> keys;
	HashMap<String, Genotype> gtByKey;
	int minSharedVariants = MIN_SHARED_VARIANTS;

	public static void main(String[] args) {
		Timer.showStdErr("Start");

		boolean debug = false;

		LikelihoodAnalysisGtPair zzz = new LikelihoodAnalysisGtPair(args);

		if (debug) {
			zzz.setDebug(debug);
			zzz.setLogLikInfoField(VCF_INFO_LOG_LIKELIHOOD);
		}

		zzz.run(debug);

		Timer.showStdErr("End");
	}

	public LikelihoodAnalysisGtPair(String args[]) {
		super(args);
		numGtAlt = 3;
		numGtNull = 2;
	}

	public LikelihoodAnalysisGtPair(String phenoCovariatesFileName, String vcfFileName) {
		super(phenoCovariatesFileName, vcfFileName);
		numGtAlt = 3;
		numGtNull = 2;
	}

	/**
	 * Calculate logistic regression's null model (or retrieve it form a cache)
	 */
	@Override
	protected double calcNullModel(int countSkip, boolean skip[], char skipChar[], double phenoNonSkip[]) {
		String skipStr = new String(skipChar);

		// Is logLikelihood cached?
		Double llNull = llNullCache.get(skipStr);
		if (llNull != null) return llNull;

		// Not found in cache? Create model and claculate it
		LogisticRegression lrNull = createNullModel(skip, countSkip, phenoNonSkip);
		this.lrNull = lrNull;
		lrNull.learn();

		llNull = lrNull.logLikelihood();
		llNullCache.put(skipStr, llNull); // Add to cache

		return llNull;
	}

	/**
	 * Keep track of the 'average' theta values (Alt model)
	 */
	protected void countModel(LogisticRegression lrAlt, LogisticRegression lrNull) {
		synchronized (thetaAltSum) {
			if (lrAlt != null) {
				double theta[] = lrAlt.getTheta();

				// Check that there are no errors
				boolean err = true;
				for (int i = 0; i < theta.length; i++)
					err &= Double.isNaN(theta[i]) || Double.isInfinite(theta[i]);

				// Add if no errors
				if (!err) {
					for (int i = 0; i < theta.length; i++)
						thetaAltSum[i] += theta[i];

					countAlt++;
				}
			}

			if (lrNull != null) {
				double theta[] = lrNull.getTheta();

				// Check that there are no errors
				boolean err = true;
				for (int i = 0; i < theta.length; i++)
					err &= Double.isNaN(theta[i]) || Double.isInfinite(theta[i]);

				// Add if no errors
				if (!err) {
					for (int i = 0; i < theta.length; i++)
						thetaNullSum[i] += theta[i];

					countNull++;
				}
			}

		}
	}

	/**
	 * Create Alt model
	 */
	LogisticRegression createAltModel(boolean skip[], int countSkip, double phenoNonSkip[], byte gti[], byte gtj[], byte gtij[]) {
		LogisticRegression lrAlt = new LogisticRegressionIrwls(numCovariates + 3); // Alt model: Include "combined" genotype gtij[]

		// Copy all covariates (except one that are skipped)
		int totalSamples = numSamples - countSkip;
		double xAlt[][] = new double[totalSamples][numCovariates + 3];

		int idx = 0;
		boolean oki = false, okj = false, okij = false;
		for (int i = 0; i < numSamples; i++) {
			if (skip[i]) continue;

			// Add genotypes
			xAlt[idx][0] = gti[i];
			xAlt[idx][1] = gtj[i];
			xAlt[idx][2] = gtij[i]; // Combined genotype: gti[i] * gtj[i]

			for (int j = 0; j < numCovariates; j++)
				xAlt[idx][j + 3] = covariates[i][j];

			if (idx > 0) {
				oki |= (xAlt[idx][0] != xAlt[idx - 1][0]);
				okj |= (xAlt[idx][1] != xAlt[idx - 1][1]);
				okij |= (xAlt[idx][2] != xAlt[idx - 1][2]);
			}

			idx++;
		}

		// Set samples
		lrAlt.setSamplesAddIntercept(xAlt, phenoNonSkip);
		lrAlt.setDebug(debug);
		setAvgThetaAltModel(lrAlt);

		this.lrAlt = lrAlt;

		if (oki && okj && okij) return lrAlt;
		return null;
	}

	/**
	 * Create null model
	 */
	LogisticRegression createNullModel(boolean skip[], int countSkip, double phenoNonSkip[], byte gti[], byte gtj[]) {
		LogisticRegression lrNull = new LogisticRegressionIrwls(numCovariates + 2); // Null model: Include "simple" genotypes

		// Copy all covariates (except one that are skipped)
		int totalSamples = numSamples - countSkip;
		double xNull[][] = new double[totalSamples][numCovariates + 2];

		int idx = 0;
		for (int i = 0; i < numSamples; i++) {
			if (skip[i]) continue;

			xNull[idx][0] = gti[i];
			xNull[idx][1] = gtj[i];

			for (int j = 0; j < numCovariates; j++)
				xNull[idx][j + 2] = covariates[i][j];

			idx++;
		}

		// Set samples
		lrNull.setSamplesAddIntercept(xNull, phenoNonSkip);
		lrNull.setDebug(debug);
		setAvgThetaNullModel(lrNull);

		this.lrNull = lrNull;
		return lrNull;
	}

	/**
	 * Are these vectors linearly dependent?
	 */
	boolean linearDependency(boolean skip[], int countSkip, byte gti[], byte gtj[], byte gtij[]) {
		int len = gti.length - countSkip;
		int n = 3;
		byte M[][] = new byte[n][len];

		// Create a matrix 'M'
		for (int i = 0, idx = 0; i < gti.length; i++) {
			if (skip[i]) continue;

			M[0][idx] = gti[i];
			M[1][idx] = gtj[i];
			M[2][idx] = gtij[i];
			idx++;
		}

		// Calculate M^t * M
		double MM[][] = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int sum = 0;

				for (int h = 0; h < len; h++)
					sum += M[i][h] * M[j][h];

				MM[i][j] = sum;
			}
		}

		// Is det( M^T * M ) zero?
		Array2DRowRealMatrix MMr = new Array2DRowRealMatrix(MM);
		double detMM = (new LUDecomposition(MMr)).getDeterminant();
		if (debug) Gpr.debug("det(MM): " + detMM + "\tMM:\n" + Gpr.toString(MM));

		return Math.abs(detMM) < EPSILON;
	}

	/**
	 * Calculate log likelihood
	 */
	public GwasResult logLikelihood(Genotype genoi, Genotype genoj) {
		String id = genoi.getId() + "-" + genoj.getId();
		byte gti[] = genoi.getGt();
		byte gtj[] = genoj.getGt();

		//---
		// Which samples should be skipped? Either missing genotype or missing phenotype
		// Also: Calculate gti * gtj vector and count number of positive entries
		//       (i.e. number of samples that have non-Ref (and non-Missing) genotypes in both variants)
		//---
		boolean skip[] = new boolean[numSamples];
		byte gtij[] = new byte[numSamples];
		int countSkip = 0, countGtij = 0;
		for (int vcfSampleNum = 0; vcfSampleNum < numSamples; vcfSampleNum++) {
			skip[vcfSampleNum] = (gti[vcfSampleNum] < 0) || (gtj[vcfSampleNum] < 0) || (pheno[vcfSampleNum] < 0);

			// Should we skip this sample?
			if (skip[vcfSampleNum]) {
				countSkip++;
			} else {
				// Calculate gt[i] * gt[j]
				gtij[vcfSampleNum] = (byte) (gti[vcfSampleNum] * gtj[vcfSampleNum]);
				if (gtij[vcfSampleNum] > 0) countGtij++; // Is it a non-Ref and non-Missing entry?
			}
		}

		//---
		// Create 'result' object
		//---
		GwasResult gwasResult = new GwasResult();
		gwasResult.genoi = genoi;
		gwasResult.genoj = genoj;
		gwasResult.gtij = gtij;

		//---
		// Sanity checks
		//---

		// No samples has both variants? Then there is not much to do.
		// To few shared varaints? We probably don't have enough statistical power anyways (not worth analysing)
		if (countGtij < minSharedVariants) {
			if (debug) Timer.show(count + "\t" + id + "\tLL_ratio: 1.0\tNot enough shared genotypes: " + countGtij);
			countModel(null);
			return gwasResult; // Not enough shared variants? Log-likelihood is probably close to zero, not worths spending time on this
		}

		// Are gti[], gtj[] and gtij[] linearly dependent?
		// If so, the model will not converge because the parameter (beta) for at least one of the gt[] will be 'NA'
		if (linearDependency(skip, countSkip, gti, gtj, gtij)) {
			if (debug) Timer.show(count + "\t" + id + "\tLL_ratio: 1.0\tLinear dependency ");
			countModel(null);
			return gwasResult; // Linear dependency? Log-likelihood is exactly zero (by definition).
		}

		//---
		// Create and fit logistic models, calculate log likelihood
		//---

		// Phenotypes without 'skipped' entries
		double phenoNonSkip[] = copyNonSkip(pheno, skip, countSkip);

		// Calculate 'Null' model (or retrieve from cache)
		LogisticRegression logRegrNull = createNullModel(skip, countSkip, phenoNonSkip, gti, gtj);
		logRegrNull.learn();
		double llNull = logRegrNull.logLikelihood();

		// Create and calculate 'Alt' model
		LogisticRegression logRegrAlt = createAltModel(skip, countSkip, phenoNonSkip, gti, gtj, gtij);
		logRegrAlt.learn();
		double llAlt = logRegrAlt.logLikelihood();

		// Calculate likelihood ratio
		double ll = 2.0 * (llAlt - llNull);

		//---
		// Save as TXT table (only used for debugging)
		//---
		if (writeToFile) {
			String idd = id.replace('/', '_');

			// ALT data
			String fileName = Gpr.HOME + "/lr_test." + idd + ".alt.txt";
			Gpr.debug("Writing 'alt data' table to :" + fileName);
			Gpr.toFile(fileName, logRegrAlt.toStringSamples());

			// NULL data
			fileName = Gpr.HOME + "/lr_test." + idd + ".null.txt";
			Gpr.debug("Writing 'null data' table to :" + fileName);
			Gpr.toFile(fileName, logRegrNull.toStringSamples());

			// ALT model
			fileName = Gpr.HOME + "/lr_test." + idd + ".alt.model.txt";
			Gpr.debug("Writing 'alt model' to :" + fileName);
			Gpr.toFile(fileName, logRegrAlt.toStringModel());
		}

		//---
		// Stats
		//---
		if (Double.isFinite(ll)) {
			// boolean show = (logLikMax < ll);
			logLikMax = Math.max(logLikMax, ll);

			if (debug) {
				// Calculate p-value
				double pval = FisherExactTest.get().chiSquareCDFComplementary(ll, deltaDf);

				Timer.show(count //
						+ "\t" + id //
						+ "\tLL_ratio: " + ll //
						+ "\tp-value: " + pval //
						+ "\tLL_alt: " + llAlt //
						+ "\tLL_null: " + llNull //
						+ "\tLL_ratio_max: " + logLikMax //
						+ (verbose ? "\n\tModel Alt  : " + logRegrAlt + "\n\tModel Null : " + logRegrNull : "") //
						);
			} else if (verbose) Timer.show(count + "\tLL_ratio: " + ll + "\t" + id);
		} else {
			// Show error
			Gpr.debug("ERROR: Likelihood ratio is infinite! ID: " + id //
					+ "\n\tLR.null : " + logRegrNull //
					+ "\n\tLR.alt  : " + logRegrAlt //
					+ "\n\tLL.null : " + llNull //
					+ "\n\tLL.alt  : " + llAlt //
					);
		}

		countModel(logRegrAlt, logRegrNull);

		// Get all data into GwasData structure
		gwasResult.logLikelihoodRatioLogReg = ll;
		gwasResult.logisticRegressionAlt = logRegrAlt;
		gwasResult.logisticRegressionNull = logRegrNull;

		return gwasResult;
	}

	@Override
	public List<VcfEntry> run(VcfFileIterator vcf, boolean createList) {
		keys = new ArrayList<>();
		gtByKey = new HashMap<String, Genotype>();

		//---
		// Read VCF: Populate list of VCF entries (single thread, used for debugging and test cases)
		//---
		int count = 1;
		for (VcfEntry ve : vcf) {
			Genotype gt = new Genotype(ve);

			// Store values
			gtByKey.put(gt.getId(), gt);
			keys.add(gt.getId());

			Gpr.showMark(count++, 1);
		}

		//---
		// Calculate likelihoods
		//---

		IntStream.range(0, keys.size()) //
		.parallel() //
		.forEach(i -> {
			for (int j = i + 1; j < keys.size(); j++) {
				String keyi = keys.get(i);
				String keyj = keys.get(j);
				Genotype gti = gtByKey.get(keyi);
				Genotype gtj = gtByKey.get(keyj);

				logLikelihood(gti, gtj);
			}
		});

		Timer.show("Done VCF file: " + gtByKey.size() + " entries");

		return new ArrayList<VcfEntry>();
	}

}
