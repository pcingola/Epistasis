package ca.mcgill.pcingola.regression;

import java.util.Arrays;

/**
 * Logistic regression
 * Model fitting by gradient descent (default)
 *
 * @author pcingola
 */
public class LogisticRegression extends Regression {

	// public static final boolean IMPLEMENTATION_CORRECT = true;

	double h[]; // Predicted outputs (model output)
	double minGradient = 0.0001;
	double eta = 1.0; // Learning (gradient)

	public LogisticRegression(int size) {
		super(size);
	}

	/**
	 * AIC: Minus twice the maximized log-likelihood plus twice the number of parameters
	 */
	public double aic() {
		return deviance() + 2 * size;
	}

	@Override
	protected double calcEnergy() {
		return -logLikelihood();
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

	/**
	 * Deviance: minus twice the maximized log-likelihood.
	 */
	public double deviance() {
		return -2.0 * logLikelihood();
	}

	public double[] getH() {
		return h;
	}

	/**
	 * Has the algorithm converged?
	 */
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
			double d = (samplesY[i] == 0 ? 1.0 - out[i] : out[i]);
			loglik += Math.log(d);
		}

		return loglik;
	}

	/**
	 * Apply model to all in[]
	 */
	@Override
	public double[] predict() {
		if (out == null) {
			out = new double[samplesX.length];
			h = new double[samplesX.length];
		}

		// Calculate model for each input
		for (int i = 0; i < numSamples; i++)
			out[i] = predict(i);

		return out;
	}

	@Override
	public double predict(double[] in) {
		double h = 0.0;

		// h = beta * in
		for (int i = 0; i < size; i++)
			h += in[i] * theta[i];

		h += theta[size]; // Last value is 'bias'

		// Calculate logit(h)
		return 1.0 / (1.0 + Math.exp(-h));
	}

	/**
	 * Predict input sample 'idx'
	 */
	protected double predict(int idx) {
		double sum = 0.0;
		double in[] = samplesX[idx];

		// h = beta * in
		for (int i = 0; i < size; i++)
			sum += in[i] * theta[i];

		sum += theta[size]; // Last value is 'bias'

		h[idx] = sum; // Store 'h'

		// Calculate logit(h)
		return 1.0 / (1.0 + Math.exp(-sum));
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
