package ca.mcgill.pcingola.regression;

import java.util.Arrays;

/**
 * Logistic regression
 * Model fitting by gradient descent (default)
 *
 * @author pcingola
 */
public class LogisticRegression extends Regression {

	double minGradient = 0.0001;
	double eta = 1.0; // Learning (gradient)
	boolean skip[]; // If set to true, samples are skipped

	public LogisticRegression(int size) {
		super(size);
	}

	@Override
	protected double calcEnergy() {
		return logLikelihood();
	}

	/**
	 * Calculate gradient
	 */
	@Override
	public double[] calcGradient() {
		predict();
		Arrays.fill(gradient, 0.0); // First guess: All parameters are zero

		int countSamples = 0;
		for (int i = 0; i < numSamples; i++) {
			if (skip != null && skip[i]) continue;

			for (int j = 0; j < dim; j++)
				gradient[j] -= (samplesY[i] - out[i]) * samplesX[i][j];

			countSamples++;
		}

		// Scale: divide by number of samples
		if (countSamples > 0) {
			for (int j = 0; j < dim; j++)
				gradient[j] /= countSamples;
		}

		return gradient;
	}

	public boolean hasConverged() {
		double sum = 0.0;

		for (int i = 0; i < theta.length; i++)
			sum += Math.abs(gradient[i]);

		return sum < minGradient;
	}

	/**
	 * Calculate log likelihood (of training data)
	 * Logarithm is in natural base ('e')
	 */
	public double logLikelihood() {
		predict();

		double loglik = 0;
		double dmin = Double.MAX_VALUE;
		for (int i = 0; i < numSamples; i++) {
			if (skip != null && skip[i]) continue;

			double d = samplesY[i] == 0 ? out[i] : 1.0 - out[i];
			dmin = Math.min(dmin, d);
			loglik += Math.log(d);
		}

		return loglik;
	}

	/**
	 * Calculate log likelihood (of training data)
	 * Logarithm is in natural base ('e')
	 */
	public double logLikelihoodNull() {
		double sum = 0;
		double h = theta[theta.length - 1];
		double o = 1.0 / (1.0 + Math.exp(-h));

		for (int i = 0; i < numSamples; i++) {
			if (skip != null && skip[i]) continue;

			sum += Math.log(samplesY[i] == 0 ? o : 1.0 - o);
		}

		return sum;
	}

	public double logLikelihoodRatio() {
		return logLikelihood() / logLikelihoodNull();
	}

	@Override
	public double predict(double[] in) {
		double h = 0.0;

		// h = beta * in
		for (int i = 0; i < size; i++)
			h += in[i] * theta[i];

		h += theta[size]; // Last value is 'bias'

		// logit(h)
		return 1.0 / (1.0 + Math.exp(-h));
	}

	@Override
	public void reset() {
		super.reset();
		if (skip != null) Arrays.fill(skip, false);
	}

	public void setEta(double eta) {
		this.eta = eta;
	}

	public void setMinGradient(double minGradient) {
		this.minGradient = minGradient;
	}

	public void setSkip(boolean[] skip) {
		this.skip = skip;
	}

	public String toStringSamples() {
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < numSamples; i++) {
			for (int j = 0; j < dim; j++)
				sb.append(samplesX[i][j] + "\t");
			sb.append(predict(samplesX[i]) + "\t");
			sb.append(samplesY[i] + "\n");
		}
		return sb.toString();
	}
}
