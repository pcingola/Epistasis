package ca.mcgill.pcingola.epistasis.gwas;

import ca.mcgill.mcb.pcingola.interval.Genome;
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

	public static double LL_SHOW_LOGREG_MODEL = 0.0; // 6.0;
	public static boolean debug = false;

	public Genotype genoi, genoj;
	public String genoiId, genojId;
	public byte gtij[]; // Genotype data used to fit the logistic regression

	public double logLikelihoodRatioLogReg = 0.0; // Log likelihood ratio from Logistic Regression model
	public double pvalueLogReg = 1.0; // P-value from log-likelihood ratio in logistic regression model

	public double logLikelihoodRatioMsa = 0.0; // Log likelihood from MSA (epistasis) model
	public double likelihoodMsaAlt = 0.0; // Likelihood from MSA (ALT model)
	public double likelihoodMsaNull = 0.0; // Likelihood from MSA (NULL model)
	public LogisticRegression logisticRegressionNull; // Logistc regression Null model
	public LogisticRegression logisticRegressionAlt; // Logistc regression Alt model

	public double bayesFactor = 0.0; // Total bayes factor
	public double log10BayesFactor = 0.0; // log10( BF )
	public double bayesFactorLogReg = 0.0; // Bayes factor for logistic regression
	public double log10BayesFactorLogReg = 0.0; // log10( BF ), only logistic regression term

	public GwasResult() {
	}

	public GwasResult(Genome genome, String line) {
		parse(genome, line);
	}

	/**
	 * Calculate "total" Bayes factor (BF_logReg * BF_MSA)
	 *
	 * @param h1 = P(theta_1 | M_1)		 This is the prior distribution (LogReg Alt model) evaluated at theta_1* (max likelihood theta_1)
	 * @param h0 = P(theta_0 | M_0)		 This is the prior distribution (LogReg Null model) evaluated at theta_0* (max likelihood theta_0)
	 */
	public double bayesFactor(double h1, double h0) {
		bayesFactor = bayesFactorLogReg(h1, h0);
		log10BayesFactorLogReg = Math.log10(bayesFactor);
		if (logLikelihoodRatioMsa > 0) bayesFactor *= Math.exp(logLikelihoodRatioMsa);
		log10BayesFactor = Math.log10(bayesFactor);
		return bayesFactor;
	}

	/**
	 * Calculate Bayes factor form Logistic regression models.
	 *
	 * Notes:
	 * 		i) The integrals are approximated using Laplace's method
	 * 		ii) We do not use a ratio of logReg.likelihoodIntegralLaplaceApproximation() method due to numerical stability issues.
	 *
	 * @param h1 = P(theta_1 | M_1)		 This is the prior distribution (Alt model) evaluated at theta_1* (max likelihood theta_1)
	 * @param h0 = P(theta_0 | M_0)		 This is the prior distribution (Null model) evaluated at theta_0* (max likelihood theta_0)
	 */
	public double bayesFactorLogReg(double h1, double h0) {
		// Calculate determinants for Hessian matrices
		double detH1 = logisticRegressionAlt.detHessian();
		double detH0 = logisticRegressionNull.detHessian();

		// Calculate (2 * pi)^((N_1 - N_0) / 2)
		int diffThetaLen = logisticRegressionAlt.getTheta().length - logisticRegressionNull.getTheta().length;
		double twopik = Math.pow(2.0 * Math.PI, diffThetaLen / 2.0);
		double diffLl = logisticRegressionAlt.logLikelihood() - logisticRegressionNull.logLikelihood();

		// Bayes factor formula
		bayesFactorLogReg = twopik * Math.sqrt(detH0 / detH1) * Math.exp(diffLl);

		return bayesFactorLogReg;
	}

	public double logLik() {
		double ll = logLikelihoodRatioLogReg;
		if (logLikelihoodRatioMsa > 0) ll += logLikelihoodRatioMsa;
		return ll;
	}

	void parse(Genome genome, String line) {
		String fields[] = line.split("\t");

		int num = 3; // First 3 fields are useless info (time, variantNumber, index)
		log10BayesFactor = parseValueDouble(fields[num++]);
		bayesFactor = Math.pow(10, log10BayesFactor);

		double logLikTotal = parseValueDouble(fields[num++]);

		log10BayesFactorLogReg = parseValueDouble(fields[num++]);
		bayesFactorLogReg = Math.pow(10, log10BayesFactorLogReg);

		pvalueLogReg = parseValueDouble(fields[num++]);
		logLikelihoodRatioLogReg = parseValueDouble(fields[num++]);

		parseValueDouble(fields[num++]); // logisticRegressionAlt.LogLik
		parseValueDouble(fields[num++]); // logisticRegressionNull.LogLik

		logLikelihoodRatioMsa = parseValueDouble(fields[num++]);
		likelihoodMsaAlt = parseValueDouble(fields[num++]);
		likelihoodMsaNull = parseValueDouble(fields[num++]);

		genoiId = fields[num++]; // Genotype_i.ID
		genoi = new Genotype(genome, genoiId);

		genojId = fields[num++]; // Genotype_j.ID
		genoj = new Genotype(genome, genojId);

		if (fields.length > num) {
			// Annotations
			genoi.setAnnotataions(fields[num++]);
			genoj.setAnnotataions(fields[num++]);

			// Logistic regression parameters
			double thetaAlt[] = parseValueDoubleArray(fields[num++]);
			logisticRegressionAlt = new LogisticRegression(thetaAlt.length - 1);
			logisticRegressionAlt.setTheta(thetaAlt);

			double thetaNull[] = parseValueDoubleArray(fields[num++]);
			logisticRegressionNull = new LogisticRegression(thetaNull.length - 1);
			logisticRegressionNull.setTheta(thetaNull);

			if (debug) Gpr.debug(genoiId + "\t" + genojId + "\t" + logLikTotal //
					+ "\n\tAlt  : " + Gpr.toString(thetaAlt) //
					+ "\n\tNull : " + Gpr.toString(thetaNull) //
			);
		}
	}

	/**
	 * Parse 'key: value' pair
	 */
	String parseValue(String field) {
		String kv[] = field.split(":");
		return kv[1].trim();
	}

	/**
	 * Parse 'key: value' as a double
	 */
	double parseValueDouble(String field) {
		return Gpr.parseDoubleSafe(parseValue(field));
	}

	/**
	 * Parse an array of doubles
	 */
	double[] parseValueDoubleArray(String field) {
		String str = parseValue(field);
		str = str.replace('[', ' ').replace(']', ' ').replace(',', '\t').trim();

		String vals[] = str.split("\t");
		double vect[] = new double[vals.length];

		for (int i = 0; i < vals.length; i++)
			vect[i] = Gpr.parseDoubleSafe(vals[i]);

		return vect;
	}

	public double pvalueLogReg() {
		int deltaDf = logisticRegressionAlt.getTheta().length - logisticRegressionNull.getTheta().length;
		pvalueLogReg = FisherExactTest.get().chiSquareCDFComplementary(logLikelihoodRatioLogReg, deltaDf);
		return pvalueLogReg;
	}

	@Override
	public String toString() {
		String additionalStr = "";

		// Show LR model
		double llt = logLik();
		if (logLikelihoodRatioLogReg >= LL_SHOW_LOGREG_MODEL) {
			additionalStr = "\t" + genoi.getAnnotataions() //
					+ "\t" + genoj.getAnnotataions() //
					+ "\tLogReg_Alt: " + Gpr.toString(logisticRegressionAlt.getTheta()) //
					+ "\tLogReg_Null: " + Gpr.toString(logisticRegressionNull.getTheta()) //
			;
		}

		if (genoi != null) genoiId = genoi.getId();
		if (genoj != null) genojId = genoj.getId();

		return "log(BF): " + log10BayesFactor //
				+ "\tll_total: " + llt //
				// Logistic regression information
				+ "\tlog(BF_LogReg): " + log10BayesFactorLogReg //
				+ "\tp-value(LogReg): " + pvalueLogReg //
				+ "\tllr_LogReg: " + logLikelihoodRatioLogReg //
				+ "\tll_LogReg_ALT: " + (logisticRegressionAlt != null ? "" + logisticRegressionAlt.logLikelihood() : "") //
				+ "\tll_LogReg_NULL: " + (logisticRegressionNull != null ? "" + logisticRegressionNull.logLikelihood() : "") //
				// Epsitasis (MSA) model
				+ "\tllr_MSA: " + logLikelihoodRatioMsa //
				+ "\tlik_MSA_ALT: " + likelihoodMsaAlt //
				+ "\tlik_MSA_NULL: " + likelihoodMsaNull //
				+ "\t" + (genoiId != null ? genoiId : "") //
				+ "\t" + (genojId != null ? genojId : "") //
				+ additionalStr //
		;

	}
}
