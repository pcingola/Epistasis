package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.StreamSupport;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.probablility.FisherExactTest;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.regression.LogisticRegression;
import ca.mcgill.pcingola.regression.LogisticRegressionIrwls;

/**
 * Logistic regression log-likelihood analysis of VCF + phenotype data
 *
 * @author pcingola
 */
public class LikelihoodAnalysis {

	public static boolean WRITE_TO_FILE = false;

	public static final int PHENO_ROW_NUMBER = 0; // Covariate number zero is phenotype

	String phenoCovariatesFileName = Gpr.HOME + "/t2d1/coEvolution/coEvolution.pheno.covariates.txt";
	String vcfFileName = Gpr.HOME + "/t2d1/vcf/eff/hm.chr1.gt.vcf";
	//	String vcfFileName = Gpr.HOME + "/t2d1/vcf/eff/z.vcf";

	boolean debug = false;
	boolean writeToFile = WRITE_TO_FILE;
	int numSamples, numCovs;
	int count = 0;
	int covariatesToNormalize[] = { 11, 12 };
	int deltaDf = 1; // Difference in degrees of freedom between Alt and Null model
	double covariates[][];
	double pheno[];
	double logLik = 0;
	double logLikMax = Double.NEGATIVE_INFINITY;
	String logLikInfoField; // If not null, an INFO field is added
	String sampleIds[];
	LogisticRegression lr;
	LogisticRegression lrAlt, lrNull;

	public static void main(String[] args) {
		Timer.showStdErr("Start");

		boolean debug = false;

		LikelihoodAnalysis zzz = new LikelihoodAnalysis(args);

		if (debug) {
			zzz.setDebug(debug);
			zzz.setLogLikInfoField("LL");
		}

		zzz.run(debug);

		Timer.showStdErr("End");
	}

	public LikelihoodAnalysis(String args[]) {
		if (args.length > 1) {
			phenoCovariatesFileName = args[0];
			vcfFileName = args[1];
		}
	}

	double[] copyNonSkip(double d[], boolean skip[], int countSkip) {
		int totalSamples = numSamples - countSkip;
		double dd[] = new double[totalSamples];

		int idx = 0;
		for (int i = 0; i < numSamples; i++)
			if (!skip[i]) dd[idx++] = d[i];

		return dd;
	}

	/**
	 * Create alternative model
	 */
	LogisticRegression createAltModel(boolean skip[], int countSkip, byte gt[], double phenoNonSkip[]) {
		LogisticRegression lrAlt = new LogisticRegressionIrwls(numCovs);

		// Copy all covariates (except one that are skipped)
		int totalSamples = numSamples - countSkip;
		double xAlt[][] = new double[totalSamples][numCovs];

		int idx = 0;
		for (int i = 0; i < numSamples; i++) {
			if (skip[i]) continue;

			// First row from phenotypes file is phenotype. But we want to predict using 'genotype', so we replace it
			xAlt[idx][0] = gt[i];

			// Copy all other covariates
			for (int j = 1; j < numCovs; j++)
				xAlt[idx][j] = covariates[i][j];

			idx++;
		}

		// Set samples
		lrAlt.setSamplesAddIntercept(xAlt, phenoNonSkip);
		lrAlt.setDebug(debug);

		return lrAlt;
	}

	/**
	 * Create null model
	 */
	LogisticRegression createNullModel(boolean skip[], int countSkip, double phenoNonSkip[]) {
		LogisticRegression lrNull = new LogisticRegressionIrwls(numCovs - 1); // Null model: No genotypes

		// Copy all covariates (except one that are skipped)
		int totalSamples = numSamples - countSkip;
		double xNull[][] = new double[totalSamples][numCovs - 1];

		int idx = 0;
		for (int i = 0; i < numSamples; i++) {
			if (skip[i]) continue;

			for (int j = 0; j < numCovs - 1; j++)
				xNull[idx][j] = covariates[i][j + 1];

			idx++;
		}

		// Set samples
		lrNull.setSamplesAddIntercept(xNull, phenoNonSkip);
		lrNull.setDebug(debug);

		return lrNull;
	}

	public double getLogLik() {
		return logLik;
	}

	public double getLogLikMax() {
		return logLikMax;
	}

	public LogisticRegression getLrAlt() {
		return lrAlt;
	}

	public LogisticRegression getLrNull() {
		return lrNull;
	}

	/**
	 * Load phenotypes and covariates
	 */
	void loadPhenoAndCovariates(String phenoCovariates) {
		Timer.showStdErr("Reading PCA data form '" + phenoCovariates + "'");

		// Read "phenotypes + covariates" file
		String lines[] = Gpr.readFile(phenoCovariates).split("\n");

		// Allocate PCA matrix
		numCovs = lines.length - 1; // Row 1 is the title, other rows are covariates
		numSamples = lines[0].split("\t").length - 1; // All, but first column are samples
		covariates = new double[numSamples][numCovs];

		// Parse
		int covariateNum = -1;
		for (String line : lines) {
			line = line.trim();
			String fields[] = line.split("\t");

			// Skip title
			if (covariateNum < 0) {
				sampleIds = new String[fields.length - 1];
				for (int i = 1; i < fields.length; i++)
					sampleIds[i - 1] = fields[i];
			} else {
				// Parse covariate values
				for (int i = 1; i < fields.length; i++)
					covariates[i - 1][covariateNum] = Gpr.parseDoubleSafe(fields[i]);
			}

			covariateNum++;
		}

		//---
		// Convert phenotype to {0,1}
		//---
		pheno = new double[numSamples];
		for (int i = 0; i < numSamples; i++)
			pheno[i] = covariates[i][PHENO_ROW_NUMBER] - 1.0;

		//---
		// Normalize covariates: sex, age
		//---
		for (int cov : covariatesToNormalize)
			normalizeCovariates(cov);
	}

	/**
	 * Fit models and calculate log likelihood test
	 */
	void logLikelihood(VcfEntry ve) {
		boolean writeToFile = this.writeToFile;

		// Which samples should be skipped?
		// i.e.: Missing genotype or missing phenotype
		// Get genotypes
		byte gt[] = ve.getGenotypesScores();
		boolean skip[] = new boolean[numSamples];
		int countSkip = 0;
		for (int vcfSampleNum = 0; vcfSampleNum < numSamples; vcfSampleNum++) {
			skip[vcfSampleNum] = (gt[vcfSampleNum] < 0) || (pheno[vcfSampleNum] < 0);
			if (skip[vcfSampleNum]) countSkip++;
		}

		// Create Null and Alt models
		double phenoNonSkip[] = copyNonSkip(pheno, skip, countSkip);
		LogisticRegression lrAlt = createAltModel(skip, countSkip, gt, phenoNonSkip);
		LogisticRegression lrNull = createNullModel(skip, countSkip, phenoNonSkip);

		//---
		// Fit logistic models
		//---
		lrNull.learn();
		lrAlt.learn();

		//---
		// Calculate likelihood ratio
		//---
		double ll = 2.0 * (lrAlt.logLikelihood() - lrNull.logLikelihood());

		if (logLikInfoField != null) ve.addInfo(logLikInfoField, "" + ll);

		// Stats
		if (Double.isFinite(ll)) {
			boolean show = (logLikMax < ll);
			logLikMax = Math.max(logLikMax, ll);

			if (show || debug) {
				// Calculate p-value
				double pval = FisherExactTest.get().chiSquareCDFComplementary(ll, deltaDf);
				Gpr.debug("TODO: Calculate and check p-value (Chi-square test): " + pval);

				writeToFile |= show;
				System.out.println(ve.toStr() //
						+ "\tLL_ratio: " + ll //
						+ "\tp-value: " + pval //
						+ "\tLL_alt: " + lrAlt.logLikelihood() //
						+ "\tLL_null: " + lrNull.logLikelihood() //
						+ "\tLL_ratio_max: " + logLikMax //
						+ "\n\tModel Alt  : " + lrAlt //
						+ "\n\tModel Null : " + lrNull //
				);
			} else Timer.show(count + "\tLL_ratio: " + ll + "\t" + ve.toStr());

		} else {
			throw new RuntimeException("Likelihood ratio is infinite!\n" + ve);
		}

		//---
		// Save as TXT table (only used for debugging)
		//---
		if (debug && writeToFile) {
			// ALT data
			String fileName = Gpr.HOME + "/lr_test." + ve.getChromosomeName() + "_" + (ve.getStart() + 1) + ".alt.txt";
			Gpr.debug("Writing 'alt data' table to :" + fileName);
			Gpr.toFile(fileName, lrAlt.toStringSamples());

			// NULL data
			fileName = Gpr.HOME + "/lr_test." + ve.getChromosomeName() + "_" + (ve.getStart() + 1) + ".null.txt";
			Gpr.debug("Writing 'null data' table to :" + fileName);
			Gpr.toFile(fileName, lrNull.toStringSamples());

			// ALT model
			fileName = Gpr.HOME + "/lr_test." + ve.getChromosomeName() + "_" + (ve.getStart() + 1) + ".alt.model.txt";
			Gpr.debug("Writing 'alt model' to :" + fileName);
			Gpr.toFile(fileName, lrAlt.toStringModel());

		}

		// Used for test cases and debugging
		this.lrNull = lrNull;
		this.lrAlt = lrAlt;
		logLik = ll;

		count++;
	}

	/**
	 * Normalize (mean 0, variance 1) a covariates' row
	 */
	void normalizeCovariates(int covNum) {
		// Calculate mean
		double sum = 0;
		for (int i = 0; i < covariates.length; i++)
			sum += covariates[i][covNum];
		double avg = sum / covariates.length;

		// Calculate stddev
		sum = 0;
		for (int i = 0; i < covariates.length; i++) {
			double d = (covariates[i][covNum] - avg);
			sum += d * d;
		}
		double stddev = Math.sqrt(sum / (covariates.length - 1));

		Timer.showStdErr("Covariate " + covNum + ", mean: " + avg + ", stddev: " + stddev);

		// Normalize
		for (int i = 0; i < covariates.length; i++)
			covariates[i][covNum] = (covariates[i][covNum] - avg) / stddev;

	}

	public void run() {
		run(false);
	}

	/**
	 * Run analysis and collect some results (for test cases)
	 */
	public List<VcfEntry> run(boolean createList) {
		//---
		// Load phenotype and covariates
		//---
		loadPhenoAndCovariates(phenoCovariatesFileName);

		//---
		// Read VCF file and perform logistic regression
		// TODO: Matrix format
		//			- Use "SnpSift allelmat" format it's much faster
		//			- We need to load the whole matrix into memory (filter out singletons / MAF < 0.1% ?)
		//			- Perform a quick check for overlapping variants between two positions before using logistic regression
		//---
		VcfFileIterator vcf = new VcfFileIterator(vcfFileName);

		//---
		// Check that sample names and sample order matches
		//---
		vcf.readHeader();
		List<String> sampleNames = vcf.getVcfHeader().getSampleNames();
		int snum = 0;
		for (String s : sampleNames) {
			if (!s.equals(sampleIds[snum])) { throw new RuntimeException("Sample names do not match:" //
					+ "\n\tSample [" + snum + "] in VCF file        :  '" + s + "'" //
					+ "\n\tSample [" + snum + "] in phenotypes file :  '" + sampleIds[snum] + "'" //
			); }
			snum++;
		}

		//---
		// Calculate for each entry in VCF file (use parallel stream
		//---
		ArrayList<VcfEntry> list = new ArrayList<VcfEntry>();
		if (createList) {
			// Process and populate list of VCF entries (single thread, used for debugging and test cases)
			for (VcfEntry ve : vcf) {
				logLikelihood(ve);
				list.add(ve);
			}
		} else StreamSupport.stream(vcf.spliterator(), true).forEach(ve -> logLikelihood(ve)); // Process (do not populate list)

		return list;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setLogLikInfoField(String logLikInfoField) {
		this.logLikInfoField = logLikInfoField;
	}

	public void setWriteToFile(boolean writeToFile) {
		this.writeToFile = writeToFile;
	}
}
