package ca.mcgill.pcingola.regression;

import java.util.Arrays;

import meshi.optimizers.SimpleStepLength;
import ca.mcgill.mcb.pcingola.util.Gpr;

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

		for (int i = 0; i < numSamples; i++)
			for (int j = 0; j < dim; j++)
				gradient[j] -= (samplesY[i] - out[i]) * samplesX[i][j];

		// Scale: divide by number of samples
		for (int j = 0; j < dim; j++)
			gradient[j] /= numSamples;

		return gradient;
	}

	public boolean hasConverged() {
		double sum = 0.0;

		for (int i = 0; i < theta.length; i++)
			sum += Math.abs(gradient[i]);

		return sum < minGradient;
	}

	/**
	 * Learn only one dimension at the time
	 */
	public void lean1D() {
		for (int i = 0; i < dim; i++)
			lean1D(i);
	}

	/**
	 * Learn only 'dimOpt' dimension 
	 */
	void lean1D(int dimOpt) {
		calcGradient();

		double sum = 0;
		for (int i = 0; i < gradient.length; i++)
			sum += Math.abs(gradient[i]);

		for (int i = 0; i < gradient.length; i++)
			if (i != dimOpt) gradient[i] = 0;
			else gradient[i] *= 0.01 / sum;

		Gpr.debug("Gradient: " + Gpr.toString(gradient) + "\tenergy: " + energy + "\t" + Gpr.toString(theta));

		try {
			SimpleStepLength lineSearch = new SimpleStepLength(this);
			lineSearch.findStepLength();
		} catch (Exception e) {
			throw new RuntimeException(e);
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
			sum += Math.log(samplesY[i] == 0 ? out[i] : 1.0 - out[i]);

		return sum;
	}

	/**
	 * Calculate log likelihood (of training data)
	 * Logarithm is in natural base ('e')
	 */
	public double logLikelihoodNull() {
		double sum = 0;
		double h = theta[theta.length - 1];
		double o = 1.0 / (1.0 + Math.exp(-h));

		for (int i = 0; i < numSamples; i++)
			sum += Math.log(samplesY[i] == 0 ? o : 1.0 - o);

		return sum;
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

	public void setEta(double eta) {
		this.eta = eta;
	}

	public void setMinGradient(double minGradient) {
		this.minGradient = minGradient;
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
