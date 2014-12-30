package ca.mcgill.pcingola.regression;

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Logistic regression
 * Model fitting by gradient descent (default)
 *
 * @author pcingola
 */
public class LogisticRegression extends Regression {

	double h[]; // Predicted outputs (model output)
	double minGradient = 0.0001;
	double loglik = Double.NaN; // Log Likelihood

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
	 * Calculate the Hessian's matrix determinant
	 */
	public double detHessian() {
		int n = theta.length;
		double H[][] = new double[n][n];

		int N = getNumSamples();
		double p[] = getOut(); // Output of logistic regression is probability

		// Calculate the Hessian matrix
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				double ss = 0;
				for (int s = 0; s < N; s++)
					ss += samplesX[s][i] * samplesX[s][j] * p[s] * (1 - p[s]);

				H[i][j] = ss;
			}

		if (debug) Gpr.debug("H:\n" + Gpr.toString(H));

		// Calculate Hessian's determinant
		Array2DRowRealMatrix Hr = new Array2DRowRealMatrix(H);
		double detH = (new LUDecomposition(Hr)).getDeterminant();
		return detH;
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
	 * Calculate an integral of the likelihood using Laplace's approximation method
	 */
	public double likelihoodIntegralLaplaceApproximation() {
		// Use Lapplace's  approximation formula
		double twopik = Math.pow(2.0 * Math.PI, theta.length / 2.0);
		double detH = detHessian();

		double intLaplace = twopik * Math.exp(logLikelihood()) * Math.sqrt(1.0 / detH);
		if (debug) Gpr.debug("Integral (Laplace): " + intLaplace + "\tdet(H): " + detH + "\tll: " + logLikelihood() + "\t(2 pi)^(K/2): " + twopik);

		return intLaplace;
	}

	/**
	 * Calculate log likelihood (of training data)
	 * Logarithm is in natural base ('e')
	 */
	public double logLikelihood() {
		if (!Double.isNaN(loglik)) return loglik; // Use cached value

		predict();

		double loglik = 0;
		for (int i = 0; i < numSamples; i++) {
			double d = (samplesY[i] == 0 ? 1.0 - out[i] : out[i]);
			loglik += Math.log(d);
		}

		return loglik;
	}

	@Override
	public void needsUpdate() {
		super.needsUpdate();
		loglik = Double.NaN;
	}

	/**
	 * Apply model to all in[]
	 */
	@Override
	public double[] predict() {
		if (!predictNeedsUpdate) return out;

		if (out == null) {
			if (samplesX == null) return null;

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
		loglik = Double.NaN;
	}

	public void setMinGradient(double minGradient) {
		this.minGradient = minGradient;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("LogisticRegression:" + super.toString());

		sb.append("\tlogLik: " + loglik);

		sb.append("\th: [");
		for (int i = 0; i < h.length; i++)
			sb.append(" " + h[i]);
		sb.append("]\t");

		return sb.toString();
	}

}
