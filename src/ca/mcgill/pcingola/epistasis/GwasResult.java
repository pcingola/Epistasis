package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * Store GWAS results here
 *
 * @author pcingola
 */
public class GwasResult {

	public String idI, idJ;
	public byte gti[], gtj[], gtij[];
	public double logLikelihoodLogReg = 0.0; // Log likelihood from Logistic Regression model
	public double logLikelihoodMsa = 0.0; // Log likelihood from MSA (epistasis) model
	public double bayesFactor = 0.0;
	public double log10BayesFactor = 0.0;
	public LogisticRegression lrNull;
	public LogisticRegression lrAlt;

	public double logLik() {
		double ll = logLikelihoodLogReg;
		if (logLikelihoodMsa > 0) ll += logLikelihoodMsa;
		return ll;
	}

	@Override
	public String toString() {
		return "ll_total: " + logLik() //
				+ "\tll_LogReg: " + logLikelihoodLogReg //
				+ "\tll_MSA: " + logLikelihoodMsa //
				+ "\t" + idI //
				+ "\t" + idJ;

	}

	/**
	 * Calculate Bayes factor (form Logistic regression models)
	 * 
	 * Notes: 
	 * 		i) The integrals are approximated using Laplace's method
	 * 		ii) We do not use a reatio of logReg.likelihoodIntegralLaplaceApproximation() method due to numerical stability issues.
	 * 
	 * @param h1 = P(theta_1 | M_1)		 This is the prior distribution (Alt model) evaluated at theta_1* (max likelihood theta_1)
	 * @param h0 = P(theta_0 | M_0)		 This is the prior distribution (Null model) evaluated at theta_0* (max likelihood theta_0)
	 */
	public double bayesFactor(double h1, double h0) {
		// Calculate determinants for Hessian matrices
		double detH1 = lrAlt.detHessian();
		double detH0 = lrNull.detHessian();

		// Calculate (2 * pi)^((N_1 - N_0) / 2)
		int diffThetaLen = lrAlt.getTheta().length - lrNull.getTheta().length;
		double twopik = Math.pow(2.0 * Math.PI, diffThetaLen / 2.0);
		double diffLl = lrAlt.logLikelihood() - lrNull.logLikelihood();

		// Bayes factor formula
		bayesFactor = twopik * Math.sqrt(detH0 / detH1) * Math.exp(diffLl);

		// Calculate log10(BF)
		log10BayesFactor = Math.log10(bayesFactor);

		// Some debugging info
		Gpr.debug("A-priory distributions not set!" //
				+ "\n\tBF             : " + bayesFactor //
				+ "\n\tlog10(BF)      : " + log10BayesFactor //
				+ "\n\tdiff.theta.len : " + diffThetaLen //
				+ "\n\tdiff.LL        : " + diffLl //
				+ "\n\tll.alt         : " + lrAlt.logLikelihood() //
				+ "\n\tll.null        : " + lrNull.logLikelihood() //
				+ "\n\tdet(H_alt)     : " + detH1 //
				+ "\n\tdet(H_null)    : " + detH0 //
		);

		return bayesFactor;
	}
}
