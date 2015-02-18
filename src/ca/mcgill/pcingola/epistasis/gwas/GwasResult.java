package ca.mcgill.pcingola.epistasis.gwas;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;

import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.probablility.FisherExactTest;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.Genotype;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * Store GWAS results here
 *
 * @author pcingola
 */
public class GwasResult {

	public static final double EPSILON = 1e-6;

	public static final int MIN_SHARED_VARIANTS = 5;
	public static double LL_SHOW_LOGREG_MODEL = 0.001; // 6.0;
	public static boolean debug = false;

	String id;
	public Genotype genoi, genoj;
	public String genoiId, genojId;
	public byte gtij[]; // Genotype data used to fit the logistic regression
	double pheno[]; // Phenotypes
	int countGtij; // Count number of samples having non-Ref and non-Missing genotypes in both variants

	int countSkip; // Number of samples skipped
	boolean skip[]; // Samples to skip (e.g. missing genotype or missing phenotype info
	String skipKey; // A string symbolizing the skipped samples. Used for caching results (null model)

	public double likelihoodLogRegAlt = 0.0; // Likelihood from logistic regression (ALT model)
	public double likelihoodLogRegNull = 0.0; // Likelihood from logistic regression (NULL model)
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

	public GwasResult(Genotype genoi, double pheno[]) {
		this.genoi = genoi;
		genoj = null;
		this.pheno = pheno;
		id = genoi.getId();
	}

	public GwasResult(Genotype genoi, Genotype genoj, double pheno[]) {
		this.genoi = genoi;
		this.genoj = genoj;
		this.pheno = pheno;
		id = genoi.getId() + "-" + genoj.getId();
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

	/**
	 * Which samples should be skipped?
	 */
	public String calcSkip() {
		if (genoj == null) return calcSkipSingle(); // Single variant

		calcSkipPair(); // Pair of variants
		return null;
	}

	/**
	 * Which samples should be skipped? Either missing genotype or missing phenotype
	 * Also: Calculate gti * gtj vector and count number of positive entries
	 * (i.e. number of samples that have non-Ref (and non-Missing) genotypes in both variants)
	 */
	void calcSkipPair() {
		// Initialize
		int numSamples = getNumSamples();
		byte gti[] = genoi.getGt();
		byte gtj[] = genoj.getGt();

		skip = new boolean[numSamples];
		gtij = new byte[numSamples];
		countSkip = 0;
		countGtij = 0;

		// Which samples should be skipped?
		gtij = new byte[numSamples];
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
	}

	/**
	 * Which samples should be skipped? Either missing genotype or missing phenotype
	 * Note: Create 'cache' key
	 */
	String calcSkipSingle() {
		int numSamples = getNumSamples();
		skip = new boolean[numSamples];
		char skipChar[] = new char[numSamples];
		countSkip = 0;
		byte[] gt = genoi.getGt();
		for (int vcfSampleNum = 0; vcfSampleNum < numSamples; vcfSampleNum++) {
			skip[vcfSampleNum] = (gt[vcfSampleNum] < 0) || (pheno[vcfSampleNum] < 0);
			if (skip[vcfSampleNum]) {
				countSkip++;
				skipChar[vcfSampleNum] = '1';
			} else skipChar[vcfSampleNum] = '0';
		}

		skipKey = new String(skipChar);
		return skipKey;
	}

	public int getCountSkip() {
		return countSkip;
	}

	public String getId() {
		return id;
	}

	int getNumSamples() {
		return genoi.numberSamples();
	}

	public boolean[] getSkip() {
		return skip;
	}

	public String getSkipKey() {
		return skipKey;
	}

	/**
	 * Are these vectors linearly dependent?
	 */
	boolean linearDependency() {
		byte gti[] = genoi.getGt();
		byte gtj[] = genoj.getGt();

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

	/**
	 * Return phenotypes (only the ones that should not be skipped)
	 */
	public double[] phenoNoSkip() {
		int totalSamples = getNumSamples() - countSkip;
		double phenoNoSkip[] = new double[totalSamples];

		int idx = 0;
		for (int i = 0; i < getNumSamples(); i++)
			if (!skip[i]) phenoNoSkip[idx++] = pheno[i];

		return phenoNoSkip;
	}

	/**
	 * Calculate p-value form logistic regression
	 */
	protected double pvalueLogReg() {
		int deltaDf = logisticRegressionAlt.getTheta().length - logisticRegressionNull.getTheta().length;
		pvalueLogReg = FisherExactTest.get().chiSquareCDFComplementary(logLikelihoodRatioLogReg, deltaDf);
		return pvalueLogReg;
	}

	/**
	 * Set logistic regression's Alt and Null models
	 */
	public void setLogRegModels(LogisticRegression logRegrAlt, LogisticRegression logRegrNull) {
		logisticRegressionAlt = logRegrAlt;
		logisticRegressionNull = logRegrNull;

		// Calculate likelihood ratio
		likelihoodLogRegAlt = logRegrAlt.logLikelihood();
		likelihoodLogRegNull = logRegrNull.logLikelihood();
		logLikelihoodRatioLogReg = 2.0 * (likelihoodLogRegAlt - likelihoodLogRegNull);

		// Calculate p-value
		pvalueLogReg();
	}

	/**
	 * Should we filter out this variant pair?
	 */
	public boolean shouldFilter() {
		// No samples has both variants? Then there is not much to do.
		// To few shared variants? We probably don't have enough statistical
		// power anyways (not worth analyzing)
		if (countGtij < MIN_SHARED_VARIANTS) {
			if (debug) Timer.show(id + "\t" + id + "\tLL_ratio: 1.0\tNot enough shared genotypes: " + countGtij);
			return true; // Not enough shared variants? Log-likelihood is probably close to zero, not worths spending time on this
		}

		// Are gti[], gtj[] and gtij[] linearly dependent?
		// If so, the model will not converge because the parameter (beta) for at least one of the gt[] will be 'NA'
		if (linearDependency()) {
			if (debug) Timer.show(id + "\tLinear dependency");
			return true; // Linear dependency? Log-likelihood is exactly zero (by definition).
		}

		// Is gti[] (or gtj[]) having non-zero entries at the same places as gtij[] ?
		// If so, the model will converge to high oposing parameters
		// between gti[] (or gtj[]) and gtij[] (probably meaningless)
		if (variantDependency()) {
			if (debug) Timer.show(id + "\tVariant match on every gtij[]");
			return true; // Linear dependency? Log-likelihood is exactly zero (by definition).
		}

		return false;
	}

	@Override
	public String toString() {
		String additionalStr = "";

		// Show LR model
		double llt = logLik();
		if (logLikelihoodRatioLogReg >= LL_SHOW_LOGREG_MODEL) {
			additionalStr = "" //
					+ "\t" + genoi.getAnnotataions() //
					+ "\t" + genoj.getAnnotataions() //
					+ "\tLogReg_Alt: " + Gpr.toString(logisticRegressionAlt.getTheta()) //
					+ "\tLogReg_Null: " + Gpr.toString(logisticRegressionNull.getTheta()) //
			;
		}

		if (genoi != null) genoiId = genoi.getId();
		if (genoj != null) genojId = genoj.getId();

		return "log10(BF): " + log10BayesFactor //
				+ "\tll_total: " + llt //
				// Logistic regression information
				+ "\tlog10(BF_LogReg): " + log10BayesFactorLogReg //
				+ "\tp-value(LogReg): " + pvalueLogReg //
				+ "\tllr_LogReg: " + logLikelihoodRatioLogReg //
				+ "\tll_LogReg_ALT: " + (logisticRegressionAlt != null ? "" + logisticRegressionAlt.logLikelihood() : "") //
				+ "\tll_LogReg_NULL: " + (logisticRegressionNull != null ? "" + logisticRegressionNull.logLikelihood() : "") //
				// Epsitasis (MSA) model
				+ "\tllr_MSA: " + logLikelihoodRatioMsa //
				+ "\tlik_MSA_ALT: " + likelihoodMsaAlt //
				+ "\tlik_MSA_NULL: " + likelihoodMsaNull //
				// Genotype IDs
				+ "\t" + (genoiId != null ? genoiId : "") //
				+ "\t" + (genojId != null ? genojId : "") //
				// Additional info
				+ additionalStr //
		;
	}

	/**
	 * Are these vectors linearly dependent?
	 */
	boolean variantDependency() {
		byte gti[] = genoi.getGt();
		byte gtj[] = genoj.getGt();

		boolean eqI = true;
		boolean eqJ = true;
		for (int i = 0; i < gti.length; i++) {
			boolean isVariantI = (gti[i] > 0);
			boolean isVariantJ = (gtj[i] > 0);
			boolean isVariantIj = (gtij[i] > 0);

			eqI &= (isVariantI == isVariantIj); // Are all variant entries in gti[] equal to variant entries in gtij[]?
			eqJ &= (isVariantJ == isVariantIj); // Are all variant entries in gtj[] equal to variant entries in gtij[]?

			// If both are different, we are done
			if (!eqI && !eqJ) return false;
		}

		// All entries in either gti[] or gtj[] are non-zero in the same positions as gtij[]
		return true;
	}
}
