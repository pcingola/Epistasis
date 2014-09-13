package ca.mcgill.pcingola.regression;

import java.util.Arrays;
import java.util.Random;

/**
 * Performs a regression
 * Single output value
 */
public abstract class Regression {

	boolean debug = false;
	int numSamples;
	int size;
	int maxIterations = 10000; // Maximum number of iterations
	double beta[]; // Parameters
	double x[][]; // Input data
	double y[]; // Output (real outputs)
	double out[]; // Predicted outputs (model output)
	boolean outputValid;
	double gradient[];
	Random rand;

	public Regression(int size) {
		this.size = size;
		beta = new double[size + 1]; // Add one parameter (e.g. intercept in linear models)
		gradient = new double[size + 1];
		initModelRand();
	}

	void checkSamples(double in[][], double out[]) {
		if (in.length != out.length) throw new RuntimeException("Sample dimensions do not match. Dim(in) = [ " + in.length + " , " + in[0].length + " ], Dim(out) = " + out.length);
		if (in[0].length != size) throw new RuntimeException("Input dimension does not model size. Dim(in) = [ " + in.length + " , " + in[0].length + " ], Dim(out) = " + size);
	}

	public double[] getBeta() {
		return beta;
	}

	public double[] getOut() {
		return out;
	}

	public double[][] getX() {
		return x;
	}

	public double[] getY() {
		return y;
	}

	/**
	 * Has this model converged?
	 */
	public abstract boolean hasConverged();

	/**
	 * Randomly initialize a model
	 */
	public void initModelRand() {
		for (int i = 0; i < size + 1; i++)
			beta[i] = randOne();

		outputValid = false;
	}

	/**
	 * Learn: Fit model
	 */
	public double[] learn() {
		outputValid = false; // Invalidate previous 'predicted' results

		int it = 0;
		do {
			if (debug) System.out.println("Iteration:  " + it);
			learnIteration();
			it++;
		} while ((it < maxIterations) && !hasConverged());

		return beta;
	}

	/**
	 * Learn parameters (single iteration)
	 */
	public abstract void learnIteration();

	/**
	 * Apply model to all in[]
	 */
	public double[] predict() {
		if (outputValid) return out;

		if (out == null) out = new double[x.length];

		// Calculate model for each input
		for (int i = 0; i < numSamples; i++)
			out[i] = predict(x[i]);

		outputValid = true;
		return out;
	}

	/**
	 * Apply model to in[]
	 */
	public abstract double predict(double in[]);

	protected double randOne() {
		double r = (rand != null ? rand.nextDouble() : Math.random());
		return 2 * r - 1;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

	public void setModel(double beta[]) {
		if (this.beta != null && this.beta.length != beta.length) throw new RuntimeException("Number of parameters does not match model size: " + this.beta.length + " != " + beta.length);
		this.beta = Arrays.copyOf(beta, beta.length);
	}

	public void setRand(Random rand) {
		this.rand = rand;
	}

	/**
	 * Data is copied into new arrays
	 * Input samples are added one column at the end with value '1' (intercept or 'bias' term in models)
	 */
	public void setSamples(double in[][], double out[]) {
		checkSamples(in, out);

		numSamples = in.length;
		x = new double[numSamples][size + 1];

		// Copy data
		for (int i = 0; i < numSamples; i++) {
			for (int j = 0; j < size; j++)
				x[i][j] = in[i][j];

			x[i][size] = 1.0; // Add 'intercept' term
		}

		y = Arrays.copyOf(out, out.length);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < size; i++) {
			if (sb.length() > 0) {
				if (beta[i] >= 0) sb.append(" + ");
				else sb.append(" ");
			}
			sb.append(beta[i] + " * in_" + i);
		}
		sb.append((beta[size] >= 0 ? " + " : " ") + beta[size]);
		return sb.toString();
	}

}
