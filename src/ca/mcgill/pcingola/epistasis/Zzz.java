package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.gwas.GwasResult;
import ca.mcgill.pcingola.epistasis.likelihood.ParameterDistributionModel;

public class Zzz {

	public static boolean debug = false;

	static double bfMsaLog10(GwasResult gr) {
		return 0;
	}

	public static void main(String[] args) {
		String file = Gpr.HOME + "/snpEff/epistasis/gwas/gwas.30.head.txt";
		String gwasTheta = Gpr.HOME + "/snpEff/epistasis/gwas/gwasCoeff_altNull.txt";

		// Load distributions
		ParameterDistributionModel pdmAlt = new ParameterDistributionModel(gwasTheta, "ALT");
		Gpr.debug(pdmAlt);

		ParameterDistributionModel pdmNull = new ParameterDistributionModel(gwasTheta, "NULL");
		Gpr.debug(pdmNull);

		// Load GwsResults
		Genome genome = new Genome("test");
		LineFileIterator lfi = new LineFileIterator(file);
		for (String line : lfi) {
			GwasResult gr = new GwasResult(genome, line);
			if (debug) Gpr.debug(gr + "\n");

			double thetaAlt[] = gr.logisticRegressionAlt.getTheta();
			double thetaNull[] = gr.logisticRegressionNull.getTheta();

			double pThetaAlt = pdmAlt.p(thetaAlt);
			double pThetaNull = pdmNull.p(thetaNull);
			double pThetaRatio = pThetaAlt / pThetaNull;
			double pThetaRatioLog10 = Math.log10(pThetaRatio);

			double bfMsaLog10 = bfMsaLog10(gr);

			// Total Bayes Factor
			double bfLog10 = gr.log10BayesFactorLogReg + pThetaRatioLog10 + bfMsaLog10;

			if (debug) {
				Gpr.debug("p(theta_alt): " + pThetaAlt //
						+ "\tp(theta_null): " + pThetaNull //
						+ "\tratio: " + pThetaRatio + "\tlog10(pThetaRatio): " + pThetaRatioLog10 //
						+ "\np_alt: " + pdmAlt.toString(thetaAlt) //
						+ "\np_null: " + pdmNull.toString(thetaNull) //
						);
			} else {
				System.out.println(bfLog10 //
						+ "\tp(theta_alt): " + pThetaAlt //
						+ "\tp(theta_null): " + pThetaNull //
						+ "\tratio: " + pThetaRatio + "\tlog10(pThetaRatio): " + Math.log10(pThetaRatio) //
						);
			}
		}
	}
}
