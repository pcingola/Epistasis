package ca.mcgill.pcingola.epistasis.likelihood;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.IntStream;

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
public class LogisticRegressionGtPair extends LogisticRegressionGt {

	ArrayList<String> keys;
	HashMap<String, Genotype> gtByKey;

	public static void main(String[] args) {
		Timer.showStdErr("Start");

		boolean debug = false;

		LogisticRegressionGtPair zzz = new LogisticRegressionGtPair(args);

		if (debug) {
			zzz.setDebug(debug);
			zzz.setLogLikInfoField(VCF_INFO_LOG_LIKELIHOOD);
		}

		zzz.init();
		zzz.run();

		Timer.showStdErr("End");
	}

	public LogisticRegressionGtPair(String args[]) {
		super(args);
		numGtAlt = 3;
		numGtNull = 2;
	}

	public LogisticRegressionGtPair(String phenoCovariatesFileName, String vcfFileName) {
		super(phenoCovariatesFileName, vcfFileName);
		numGtAlt = 3;
		numGtNull = 2;
	}

	/**
	 * Create Alt model
	 */
	@Override
	protected LogisticRegression createAltModel(GwasResult gwasResult, double phenoNonSkip[]) {
		LogisticRegression lrAlt = new LogisticRegressionIrwls(numCovariates + 3); // Alt model: Include "combined" genotype gtij[]

		// Copy all covariates (except one that are skipped)
		int totalSamples = numSamples - gwasResult.getCountSkip();
		double xAlt[][] = new double[totalSamples][numCovariates + 3];

		int idx = 0;
		boolean oki = false, okj = false, okij = false;
		byte gti[] = gwasResult.genoi.getGt();
		byte gtj[] = gwasResult.genoj.getGt();
		byte gtij[] = gwasResult.gtij;
		boolean skip[] = gwasResult.getSkip();

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

		this.lrAlt = lrAlt;

		if (oki && okj && okij) return lrAlt;
		return null;
	}

	/**
	 * Create null model
	 */
	protected LogisticRegression createNullModel(GwasResult gwasResult, double phenoNonSkip[]) {
		LogisticRegression lrNull = new LogisticRegressionIrwls(numCovariates + 2); // Null model: Include "simple" genotypes

		// Copy all covariates (except one that are skipped)
		int totalSamples = numSamples - gwasResult.getCountSkip();
		double xNull[][] = new double[totalSamples][numCovariates + 2];

		int idx = 0;
		byte gti[] = gwasResult.genoi.getGt();
		byte gtj[] = gwasResult.genoj.getGt();
		boolean skip[] = gwasResult.getSkip();

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

		this.lrNull = lrNull;
		return lrNull;
	}

	/**
	 * Calculate log likelihood
	 */
	public GwasResult logLikelihood(Genotype genoi, Genotype genoj) {
		//---
		// Create 'result' object
		//---
		GwasResult gwasResult = new GwasResult(genoi, genoj, pheno);

		// Should we filter this pair out?
		gwasResult.calcSkip();
		if (gwasResult.shouldFilter()) return gwasResult;
		Gpr.debug("CALC LOG REG!");

		//---
		// Create and fit logistic models, calculate log likelihood
		//---

		// Phenotypes without 'skipped' entries
		double phenoNonSkip[] = gwasResult.phenoNoSkip();

		// Calculate 'Null' model (or retrieve from cache)
		//LogisticRegression logRegrNull = createNullModel(skip, countSkip, phenoNonSkip, gti, gtj);
		LogisticRegression logRegrNull = createNullModel(gwasResult, phenoNonSkip);
		logRegrNull.learn();
		double llNull = logRegrNull.logLikelihood();

		// Create and calculate 'Alt' model
		LogisticRegression logRegrAlt = createAltModel(gwasResult, phenoNonSkip);
		logRegrAlt.learn();
		double llAlt = logRegrAlt.logLikelihood();

		// Calculate likelihood ratio
		double ll = 2.0 * (llAlt - llNull);

		//---
		// Save as TXT table (only used for debugging)
		//---
		if (writeToFile) {
			String idd = gwasResult.getId().replace('/', '_');

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
			logLikMax = Math.max(logLikMax, ll);

			if (debug) {
				// Calculate p-value
				double pval = FisherExactTest.get().chiSquareCDFComplementary(ll, deltaDf);

				Timer.show(count //
						+ "\t" + gwasResult.getId() //
						+ "\tLL_ratio: " + ll //
						+ "\tp-value: " + pval //
						+ "\tLL_alt: " + llAlt //
						+ "\tLL_null: " + llNull //
						+ "\tLL_ratio_max: " + logLikMax //
						+ (verbose ? "\n\tModel Alt  : " + logRegrAlt + "\n\tModel Null : " + logRegrNull : "") //
				);
			} else if (verbose) Timer.show(count + "\tLL_ratio: " + ll + "\t" + gwasResult.getId());
		} else {
			// Logitic regression is infinite: Show error
			Gpr.debug("ERROR: Likelihood ratio is infinite! ID: " + gwasResult.getId() //
					+ "\n\tLR.null : " + logRegrNull //
					+ "\n\tLR.alt  : " + logRegrAlt //
					+ "\n\tLL.null : " + llNull //
					+ "\n\tLL.alt  : " + llAlt //
			);
		}

		// Get all data into GwasData structure
		gwasResult.logLikelihoodRatioLogReg = ll;
		gwasResult.logisticRegressionAlt = logRegrAlt;
		gwasResult.logisticRegressionNull = logRegrNull;

		return gwasResult;
	}

	@Override
	public void run() {
		//---
		// Read VCF file
		//---
		List<Genotype> gts = new ArrayList<Genotype>(); // Store genotypes for split_i
		Timer.showStdErr("Reading vcf file '" + vcfFileName + "'");
		VcfFileIterator vcf = new VcfFileIterator(vcfFileName);
		for (VcfEntry ve : vcf) {
			// Store VCF entry
			Genotype geno = new Genotype(ve);
			gts.add(geno);
		}

		//---
		// Calculate likelihoods
		//---
		int count = 1;
		for (int idxi = 0; idxi < gts.size(); idxi++) {
			for (int idxj = idxi + 1; idxj < gts.size(); idxj++) {
				Genotype gti = gts.get(idxi);
				Genotype gtj = gts.get(idxj);
				if (verbose) System.out.println(gti + "\t" + gtj);

				GwasResult gwasResult = logLikelihood(gti, gtj);
				if (verbose) System.out.println(gwasResult);
				else Gpr.showMark(count++, 100);
			}
		}
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
