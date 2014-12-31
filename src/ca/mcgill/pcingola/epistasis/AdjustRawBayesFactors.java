package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.gwas.GwasResult;
import ca.mcgill.pcingola.epistasis.likelihood.ParameterDistributionModel;

/**
 * Adjust raw Bayes factors using priors distributions
 * and Ln[ Odds[ LL(MSA) ] ] empirical fomrmula
 *
 * @author pcingola
 */
public class AdjustRawBayesFactors {

	public static boolean debug = false;
	public static boolean quiet = true;
	public static final double LOG_10 = Math.log(10);

	// This models and the parameters (alpha and beta) are found empirically
	// We used real pairs of proteins "in contact" from PDB to approximate
	// the Ln[Odds] formula
	public static final double LL_MSA_APPROX_ALPHA = 0.1947;
	public static final double LL_MSA_APPROX_BETA = -1.0174;
	public static final double LL_MSA_MAX_LOG_ODDS = 6.0; // Cap LL_MSA to this value (largest seen in real data is arround 3) 

	String gwasThetaDistributionFile;
	String gwasResutsFile;

	public AdjustRawBayesFactors(String gwasThetaDistributionFile, String gwasResutsFile) {
		this.gwasThetaDistributionFile = gwasThetaDistributionFile;
		this.gwasResutsFile = gwasResutsFile;
	}

	double bfMsaLog10(GwasResult gr) {
		double llRatio = gr.logLikelihoodRatioMsa;
		double logOdds = Math.exp(LL_MSA_APPROX_ALPHA * llRatio + LL_MSA_APPROX_BETA);
		logOdds /= LOG_10; // Convert from Ln to Log10
		if (logOdds > 10) Gpr.debug("logOdds: " + logOdds + "\t" + gr.logLikelihoodRatioMsa);

		logOdds = Math.min(logOdds, LL_MSA_MAX_LOG_ODDS);
		return logOdds;
	}

	public void run() {
		//---
		// Load distributions
		//---
		ParameterDistributionModel pdmAlt = new ParameterDistributionModel(gwasThetaDistributionFile, "ALT");
		Gpr.debug(pdmAlt);

		ParameterDistributionModel pdmNull = new ParameterDistributionModel(gwasThetaDistributionFile, "NULL");
		Gpr.debug(pdmNull);

		//---
		// Read gwasResutls and adjust bayes factors
		//---
		Genome genome = new Genome("test");
		LineFileIterator lfi = new LineFileIterator(gwasResutsFile);
		for (String line : lfi) {
			GwasResult gr = new GwasResult(genome, line);
			if (debug) Gpr.debug(gr + "\n");

			double thetaAlt[] = gr.logisticRegressionAlt.getTheta();
			double thetaNull[] = gr.logisticRegressionNull.getTheta();

			double pThetaAlt = pdmAlt.p(thetaAlt);
			double pThetaNull = pdmNull.p(thetaNull);
			double pThetaRatio = pThetaAlt / pThetaNull;
			double pThetaRatioLog10 = Math.log10(pThetaRatio);

			// Update logistic regression value
			gr.log10BayesFactorLogReg += pThetaRatioLog10;

			// Update LogLikMsa (using empirical formula)
			double bfMsaLog10 = bfMsaLog10(gr);
			gr.logLikelihoodRatioMsa = bfMsaLog10;

			// Update total Bayes Factor
			double bfLog10 = gr.log10BayesFactorLogReg + bfMsaLog10;
			gr.log10BayesFactor = bfLog10;
			gr.bayesFactor = Math.pow(10, bfLog10);

			if (debug) {
				Gpr.debug("p(theta_alt): " + pThetaAlt //
						+ "\tp(theta_null): " + pThetaNull //
						+ "\tratio: " + pThetaRatio + "\tlog10(pThetaRatio): " + pThetaRatioLog10 //
						+ "\np_alt: " + pdmAlt.toString(thetaAlt) //
						+ "\np_null: " + pdmNull.toString(thetaNull) //
				);
			} else if (!quiet) {
				System.out.println(bfLog10 //
						+ "\tbfMsaLog10: " + bfMsaLog10 //
						+ "\tlog10(pThetaRatio): " + pThetaRatioLog10 //
						+ "\t" + gr //
				);
			}
		}
	}
}
