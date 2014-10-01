package ca.mcgill.pcingola.regression;

import java.util.Arrays;

/**
 * Logistic regression
 * Model fitting by gradient descent (default)
 *
 * @author pcingola
 */
public class LogisticRegression extends Regression {

	public static final boolean IMPLEMENTATION_CORRECT = false;

	double minGradient = 0.0001;
	double eta = 1.0; // Learning (gradient)

	public LogisticRegression(int size) {
		super(size);
	}

	@Override
	protected double calcEnergy() {
		return (IMPLEMENTATION_CORRECT ? -1 : 1) * (1.0 / numSamples) * logLikelihood();
	}

	/**
	 * Calculate gradient
	 * Note: Since we want the maximum likelihood and we have a minimization
	 *       algorithm, we then minimize "negative log-likelihood".
	 *       That's why the gradient is negative.
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
		for (int i = 0; i < numSamples; i++) {
			if (skip != null && skip[i]) continue;

			double d;
			if (IMPLEMENTATION_CORRECT) d = (samplesY[i] == 0 ? 1.0 - out[i] : out[i]);
			else d = (samplesY[i] == 0 ? out[i] : 1.0 - out[i]);

			loglik += Math.log(d);
		}

		return loglik;
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

}
