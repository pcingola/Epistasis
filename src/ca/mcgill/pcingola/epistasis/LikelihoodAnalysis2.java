package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.regression.LogisticRegression;
import ca.mcgill.pcingola.regression.LogisticRegressionIrwls;

/**
 * Logistic regression log-likelihood analysis of 2 VCF entries + phenotype data
 *
 * @author pcingola
 */
public class LikelihoodAnalysis2 extends LikelihoodAnalysis {

	public static final int MIN_SHARED_VARIANTS = 5;

	ArrayList<String> keys;
	HashMap<String, byte[]> gtByKey;
	int minSharedVariants = MIN_SHARED_VARIANTS;

	public static void main(String[] args) {
		Timer.showStdErr("Start");

		boolean debug = false;

		LikelihoodAnalysis2 zzz = new LikelihoodAnalysis2(args);

		if (debug) {
			zzz.setDebug(debug);
			zzz.setLogLikInfoField(VCF_INFO_LOG_LIKELIHOOD);
		}

		zzz.run(debug);

		Timer.showStdErr("End");
	}

	public LikelihoodAnalysis2(String args[]) {
		super(args);
	}

	protected double logLikelihood(String keyi, byte gti[], String keyj, byte gtj[]) {

		//---
		// Which samples should be skipped? Either missing genotype or missing phenotype
		// Also: Calculate gti * gtj vector and count number of positive entries (i.e. number of samples that have non-Ref (and non-Missing) genotypes in both variants)
		//---
		boolean skip[] = new boolean[numSamples];
		byte gtij[] = new byte[numSamples];
		char skipChar[] = new char[numSamples];
		int countSkip = 0, countGtij = 0;
		for (int vcfSampleNum = 0; vcfSampleNum < numSamples; vcfSampleNum++) {
			skip[vcfSampleNum] = (gti[vcfSampleNum] < 0) || (gtj[vcfSampleNum] < 0) || (pheno[vcfSampleNum] < 0);

			// Should we skip this sample?
			if (skip[vcfSampleNum]) {
				countSkip++;
				skipChar[vcfSampleNum] = '1';
			} else {
				// Calculate gt[i] * gt[j]
				gtij[vcfSampleNum] = (byte) (gti[vcfSampleNum] * gtj[vcfSampleNum]);
				if (gtij[vcfSampleNum] > 0) countGtij++; // Is it a non-Ref and non-Missing entry? 

				skipChar[vcfSampleNum] = '0';
			}
		}

		// No samples has both variants? Then there is not much to do.
		// To few shared varaints? We probably don't have enough statistical power anyways (not worth analysing)
		if (countGtij < minSharedVariants) return 0.0;

		//---
		// Create models and calculate likelihoods
		//---

		return 1.0;
	}

	/**
	 * Create null model
	 */
	@Override
	LogisticRegression createNullModel(boolean skip[], int countSkip, double phenoNonSkip[]) {
		LogisticRegression lrNull = new LogisticRegressionIrwls(numCovariates - 1); // Null model: No genotypes

		// Copy all covariates (except one that are skipped)
		int totalSamples = numSamples - countSkip;
		double xNull[][] = new double[totalSamples][numCovariates - 1];

		int idx = 0;
		for (int i = 0; i < numSamples; i++) {
			if (skip[i]) continue;

			for (int j = 0; j < numCovariates - 1; j++)
				xNull[idx][j] = covariates[i][j + 1];

			idx++;
		}

		// Set samples
		lrNull.setSamplesAddIntercept(xNull, phenoNonSkip);
		lrNull.setDebug(debug);

		this.lrNull = lrNull;
		return lrNull;
	}

	@Override
	public List<VcfEntry> run(VcfFileIterator vcf, boolean createList) {
		keys = new ArrayList<>();
		gtByKey = new HashMap<String, byte[]>();

		//---
		// Read VCF: Populate list of VCF entries (single thread, used for debugging and test cases)
		//---
		Timer.show("Reading VCF file");
		int count = 1;
		for (VcfEntry ve : vcf) {
			String key = ve.toStr();
			byte gt[] = ve.getGenotypesScores();

			// Store values
			gtByKey.put(key, gt);
			keys.add(key);

			Gpr.showMark(count++, 1000);
		}

		//---
		// Calculate likelihoods
		//---

		for (int i = 0; i < keys.size(); i++) {
			for (int j = i + 1; j < keys.size(); j++) {
				String keyi = keys.get(i);
				String keyj = keys.get(j);
				byte gti[] = gtByKey.get(keyi);
				byte gtj[] = gtByKey.get(keyj);

				logLikelihood(keyi, gti, keyj, gtj);
			}
		}

		Timer.show("Done VCF file: " + gtByKey.size() + " entries");

		return new ArrayList<VcfEntry>();
	}
}
