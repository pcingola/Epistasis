package ca.mcgill.pcingola.epistasis.gwas;

import ca.mcgill.mcb.pcingola.probablility.FisherExactTest;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.Genotype;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * Store GWAS results here
 *
 * @author pcingola
 */
public class GwasResult {

	public static double LL_SHOW_LOGREG_MODEL = 6.0;

	public Genotype genoi, genoj;
	public byte gtij[]; // Genotype data used to fit the logistic regression

	public double logLikelihoodLogReg = 0.0; // Log likelihood from Logistic Regression model
	public double pvalueLogReg = 1.0; // P-value from log-likelihood in logistic regression model

	public double logLikelihoodMsa = 0.0; // Log likelihood from MSA (epistasis) model
	public double likelihoodMsaAlt = 0.0; // Likelihood from MSA (ALT model)
	public double likelihoodMsaNull = 0.0; // Likelihood from MSA (NULL model)
	public LogisticRegression lrNull; // Logistc regression Null model
	public LogisticRegression lrAlt; // Logistc regression Alt model

	public double bayesFactorLogReg = 0.0; // Bayes factor for logistic regression
	public double bayesFactor = 0.0; // Total bayes factor
	public double log10BayesFactor = 0.0; // log10( BF )

	/**
	 * Calculate "total" Bayes factor (BF_logReg * BF_MSA)
	 *
	 * @param h1 = P(theta_1 | M_1)		 This is the prior distribution (LogReg Alt model) evaluated at theta_1* (max likelihood theta_1)
	 * @param h0 = P(theta_0 | M_0)		 This is the prior distribution (LogReg Null model) evaluated at theta_0* (max likelihood theta_0)
	 */
	public double bayesFactor(double h1, double h0) {
		bayesFactor = bayesFactorLogReg(h1, h0);
		if (logLikelihoodMsa > 0) bayesFactor *= Math.exp(logLikelihoodMsa);
		log10BayesFactor = Math.log10(bayesFactor);
		return bayesFactor;
	}

	/**
	 * Calculate Bayes factor form Logistic regression models.
	 *
	 * Notes:
	 * 		i) The integrals are approximated using Laplace's method
	 * 		ii) We do not use a reatio of logReg.likelihoodIntegralLaplaceApproximation() method due to numerical stability issues.
	 *
	 * @param h1 = P(theta_1 | M_1)		 This is the prior distribution (Alt model) evaluated at theta_1* (max likelihood theta_1)
	 * @param h0 = P(theta_0 | M_0)		 This is the prior distribution (Null model) evaluated at theta_0* (max likelihood theta_0)
	 */
	public double bayesFactorLogReg(double h1, double h0) {
		// Calculate determinants for Hessian matrices
		double detH1 = lrAlt.detHessian();
		double detH0 = lrNull.detHessian();

		// Calculate (2 * pi)^((N_1 - N_0) / 2)
		int diffThetaLen = lrAlt.getTheta().length - lrNull.getTheta().length;
		double twopik = Math.pow(2.0 * Math.PI, diffThetaLen / 2.0);
		double diffLl = lrAlt.logLikelihood() - lrNull.logLikelihood();

		// Bayes factor formula
		bayesFactorLogReg = twopik * Math.sqrt(detH0 / detH1) * Math.exp(diffLl);

		return bayesFactorLogReg;
	}

	public double logLik() {
		double ll = logLikelihoodLogReg;
		if (logLikelihoodMsa > 0) ll += logLikelihoodMsa;
		return ll;
	}

	public double pvalueLogReg() {
		int deltaDf = lrAlt.getTheta().length - lrNull.getTheta().length;
		pvalueLogReg = FisherExactTest.get().chiSquareCDFComplementary(logLikelihoodLogReg, deltaDf);
		return pvalueLogReg;
	}

	@Override
	public String toString() {
		String additionalStr = "";

		// Show LR model
		double llt = logLik();
		if (logLikelihoodLogReg >= LL_SHOW_LOGREG_MODEL) {
			additionalStr = "\t" + genoi.getAnnotataions() //
					+ "\t" + genoj.getAnnotataions() //
					+ "\n\tAlt  : " + Gpr.toString(lrAlt.getTheta()) //
					+ "\n\tNull :             " + Gpr.toString(lrNull.getTheta()) //

			;
		}

		return "log(BF): " + log10BayesFactor //
				+ "\tp-value(LogReg): " + pvalueLogReg //
				+ "\tll_total: " + llt //
				+ "\tll_LogReg: " + logLikelihoodLogReg //
				+ "\tll_MSA: " + logLikelihoodMsa //
				+ "\tlik_MSA_ALT: " + likelihoodMsaAlt //
				+ "\tlik_MSA_NULL: " + likelihoodMsaNull //
				+ "\t" + genoi.getId() //
				+ "\t" + genoj.getId() //
				+ additionalStr //
		;

	}
}
