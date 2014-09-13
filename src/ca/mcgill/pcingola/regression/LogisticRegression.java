package ca.mcgill.pcingola.regression;

import java.util.Arrays;

/**
 * Logistic regression
 * Model fitting by gradient descent
 *
 * @author pcingola
 */
public class LogisticRegression extends Regression {

	double minGradient = 0.0001;
	double eta = 1.0; // Learning (gradient)

	public LogisticRegression(int size) {
		super(size);
	}

	/**
	 * Calculate gradient
	 */
	public double[] gradient() {
		predict();
		Arrays.fill(gradient, 0.0);

		int dim = size + 1;
		for (int i = 0; i < numSamples; i++)
			for (int j = 0; j < dim; j++)
				gradient[j] += (y[i] - out[i]) * x[i][j];

		// Scale: divide by number of samples
		for (int j = 0; j < dim; j++)
			gradient[j] /= numSamples;

		return gradient;
	}

	@Override
	public boolean hasConverged() {
		double sum = 0.0;
		for (int i = 0; i < beta.length; i++)
			sum += Math.abs(gradient[i]);
		return sum < minGradient;
	}

	@Override
	public void learnIteration() {
		outputValid = false;
		predict();
		gradient();

		for (int i = 0; i < beta.length; i++) {
			if (debug) System.out.println("\t" + i + "\tbeta: " + beta[i] + "\tgradient: " + gradient[i]);
			beta[i] += gradient[i];
		}

	}

	/**
	 * Calculate log likelihood (of training data)
	 * Logarithm is in natural base ('e')
	 */
	public double logLikelihood() {
		predict();

		double sum = 0;
		for (int i = 0; i < numSamples; i++)
			sum += Math.log(y[i] == 0 ? out[i] : 1.0 - out[i]);

		return sum;
	}

	/**
	 * Calculate log likelihood (of training data)
	 * Logarithm is in natural base ('e')
	 */
	public double logLikelihoodNull() {
		double sum = 0;
		double h = beta[beta.length - 1];
		double o = 1.0 / (1.0 + Math.exp(-h));

		for (int i = 0; i < numSamples; i++)
			sum += Math.log(y[i] == 0 ? o : 1.0 - o);

		return sum;
	}

	@Override
	public double predict(double[] in) {
		double h = 0.0;

		// beta * in
		for (int i = 0; i < size; i++)
			h += in[i] * beta[i];

		h += beta[size]; // Last value is 'bias'

		return 1.0 / (1.0 + Math.exp(-h));
	}

	public void setEta(double eta) {
		this.eta = eta;
	}

	public void setMinGradient(double minGradient) {
		this.minGradient = minGradient;
	}

}
