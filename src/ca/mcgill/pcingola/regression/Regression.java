package ca.mcgill.pcingola.regression;

import java.util.Arrays;
import java.util.Random;

import meshi.optimizers.Energy;
import meshi.optimizers.Minimizer;
import meshi.optimizers.SteepestDecent;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Generic regression for single values output
 */
public abstract class Regression extends Energy {

	boolean debug = false;
	boolean skip[]; // If set to true, samples are skipped
	int numSamples;
	int size;
	int maxIterations = 10000; // Maximum number of iterations
	double samplesX[][]; // Samples: Input data
	double samplesY[]; // Samples: Output (real outputs)
	double out[]; // Predicted outputs (model output)
	Random rand;
	Minimizer minimizer;

	public Regression(int size) {
		super(size + 1); // Add one for intercept
		this.size = size;
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

	public boolean[] getSkip() {
		return skip;
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
		if (rand != null) initModelRand();
		else Arrays.fill(theta, 0.0);

		// if (minimizer == null) minimizer = new GradientDecent(this);
		if (minimizer == null) minimizer = new SteepestDecent(this);
		minimizer.setDebug(debug);
		minimizer.run();
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
		double r = rand != null ? rand.nextDouble() : Math.random();
		return 2 * r - 1;
	}

	@Override
	public void reset() {
		super.reset();
		minimizer = null;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

	public void setMinnimizer(Minimizer minnimizer) {
		minimizer = minnimizer;
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

	public void setSkip(boolean[] skip) {
		this.skip = skip;
	}

	@Override
	public String toString() {
		return energy + ", model: " + Gpr.toString(theta);
	}

	public String toStringSamples() {
		StringBuilder sb = new StringBuilder();

		// Title line
		for (int j = 0; j < dim; j++)
			sb.append("in" + j + "\t");
		sb.append("predict\tout\n");

		// Show data
		for (int i = 0; i < numSamples; i++) {
			if (skip[i]) continue; // Don't show skipped samples

			for (int j = 0; j < dim; j++)
				sb.append(samplesX[i][j] + "\t");

			sb.append(predict(samplesX[i]) + "\t");
			sb.append(samplesY[i] + "\n");
		}

		return sb.toString();
	}

}
