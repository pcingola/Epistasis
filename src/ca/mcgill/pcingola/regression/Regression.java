package ca.mcgill.pcingola.regression;

import java.util.Random;

import meshi.optimizers.Energy;
import meshi.optimizers.Minimizer;
import meshi.optimizers.SteepestDecent;

/**
 * Generic regression for single values output
 */
public abstract class Regression extends Energy {

	boolean debug = false;
	int numSamples;
	int size;
	int maxIterations = 10000; // Maximum number of iterations
	double samplesX[][]; // Samples: Input data
	double samplesY[]; // Samples: Output (real outputs)
	double out[]; // Predicted outputs (model output)
	Random rand;
	Minimizer minnimizer;

	public Regression(int size) {
		super(size + 1); // Add one for intercept
		this.size = size;
		initModelRand();
	}

	void checkSamples(double in[][], double out[]) {
		if (in.length != out.length) throw new RuntimeException("Sample dimensions do not match. Dim(in) = [ " + in.length + " , " + in[0].length + " ], Dim(out) = " + out.length);
		if (in[0].length != (dim - 1)) throw new RuntimeException("Input dimension does not model size. Dim(in) = [ " + in.length + " , " + in[0].length + " ], Dim(out) = " + (dim - 1));
	}

	public double[] getOut() {
		return out;
	}

	public double[][] getSamplesX() {
		return samplesX;
	}

	public double[] getSamplesY() {
		return samplesY;
	}

	/**
	 * Randomly initialize a model
	 */
	public void initModelRand() {
		for (int i = 0; i < dim; i++)
			theta[i] = randOne();

		needsUpdate();
	}

	/**
	 * Learn: Fit model
	 */
	public double[] learn() {
		if (minnimizer == null) minnimizer = new SteepestDecent(this);
		minnimizer.run();
		return theta;
	}

	/**
	 * Apply model to all in[]
	 */
	public double[] predict() {
		if (out == null) out = new double[samplesX.length];

		// Calculate model for each input
		for (int i = 0; i < numSamples; i++)
			out[i] = predict(samplesX[i]);

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

	public void setMinnimizer(Minimizer minnimizer) {
		this.minnimizer = minnimizer;
	}

	public void setModel(double theta[]) {
		if (this.theta != null && this.theta.length != theta.length) throw new RuntimeException("Number of parameters does not match model size: " + this.theta.length + " != " + theta.length);
		System.arraycopy(theta, 0, this.theta, 0, theta.length);
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
		samplesX = new double[numSamples][dim + 1];

		// Copy data
		for (int i = 0; i < numSamples; i++) {
			for (int j = 0; j < size; j++)
				samplesX[i][j] = in[i][j];

			samplesX[i][size] = 1.0; // Add 'intercept' term
		}

		if (samplesY == null) samplesY = new double[out.length];
		System.arraycopy(out, 0, samplesY, 0, out.length);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(energy + ", model: ");
		for (int i = 0; i < size; i++) {
			if (sb.length() > 0) {
				if (theta[i] > 0 && i > 0) sb.append(" + ");
				else sb.append(" ");
			}
			sb.append(theta[i] + " * in_" + i);
		}
		sb.append((theta[size] >= 0 ? " + " : " ") + theta[size]);
		return sb.toString();
	}
}
