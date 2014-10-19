package meshi.optimizers;

import java.util.Arrays;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A generic class for all energy functions to be minimized.
 */
public abstract class Energy {

	protected int dim;
	protected double[] theta; // Parameters to optimize
	protected double[] thetaBest; // Best parameters so far
	protected double[] gradient; // Parameter's gradient
	protected double energy;
	protected boolean energyNeedsUpdate = true;
	protected boolean gradientNeedsUpdate = true;

	public Energy(int dim) {
		this.dim = dim;
		theta = new double[dim];
		thetaBest = new double[dim];
		gradient = new double[dim];
		energy = Double.NaN;
	}

	public void addTheta(int idx, double addToXi) {
		theta[idx] += addToXi;
		needsUpdate();
	}

	public void addThetaBestGradient(double alpha) {
		for (int i = 0; i < theta.length; i++)
			theta[i] = thetaBest[i] + alpha * gradient[i];

		needsUpdate();
	}

	public void addThetaGradient(double alpha) {
		for (int i = 0; i < theta.length; i++)
			theta[i] += alpha * gradient[i];

		needsUpdate();
	}

	/**
	 * Calculate energy
	 */
	protected abstract double calcEnergy();

	/**
	 * Calculate gradient
	 */
	protected abstract double[] calcGradient();

	/**
	 * Copy current gradient[] to an array
	 */
	public void copyGradient(double copyGrad[]) {
		System.arraycopy(gradient, 0, copyGrad, 0, gradient.length);
	}

	/**
	 * Copy current x[] to an array
	 */
	public void copyTheta(double copyX[]) {
		System.arraycopy(theta, 0, copyX, 0, theta.length);
	}

	/**
	 * Evaluate energy (store value for future use)
	 */
	public double evaluate() {
		updateEnergy();
		updateGradient();
		return energy;
	}

	/**
	 * Returns latest energy value
	 **/
	public double getEnergy() {
		return energy;
	}

	/**
	 * Get latest gradient
	 */
	public double[] getGradient() {
		return gradient;
	}

	/**
	 * Scaled gradient
	 */
	public double[] getGradient(double factor) {
		double ng[] = new double[gradient.length];
		for (int i = 0; i < ng.length; i++)
			ng[i] = factor * gradient[i];
		return ng;
	}

	/**
	 * Scaled gradient
	 */
	public double[] getGradient(double factor, double buffer[]) {
		for (int i = 0; i < gradient.length; i++)
			buffer[i] = factor * gradient[i];
		return buffer;
	}

	public double[] getTheta() {
		return theta;
	}

	public double[] getThetaBest() {
		return thetaBest;
	}

	/**
	 * Gradient's norm
	 */
	public double gradientNorm() {
		double sum = 0;

		for (int i = 0; i < theta.length; i++) {
			double g = gradient[i];
			sum += g * g;
		}

		return Math.sqrt(sum);
	}

	public void needsUpdate() {
		energyNeedsUpdate = gradientNeedsUpdate = true;
	}

	public void reset() {
		needsUpdate();
		Arrays.fill(theta, 0.0);
		Arrays.fill(thetaBest, 0.0);
		Arrays.fill(gradient, 0.0);
	}

	public void setTheta(double newX[]) {
		System.arraycopy(newX, 0, theta, 0, theta.length);
		needsUpdate();
	}

	public void setTheta(int idx, double xi) {
		theta[idx] = xi;
		needsUpdate();
	}

	/**
	 * Best energy so far? Keep a copy
	 */
	public void setThetaBest() {
		System.arraycopy(theta, 0, thetaBest, 0, theta.length);
	}

	public int size() {
		return theta.length;
	}

	@Override
	public String toString() {
		return "Energy: " + (energyNeedsUpdate ? "[Needs update]" : "") //
				+ energy //
				//
				+ "\ttheta: " + Gpr.toString(theta) //
				//
				+ "\tgradient[theta] " //
				+ (gradientNeedsUpdate ? "[Needs update]" : "") //
				+ ": " + Gpr.toString(gradient);
	}

	/**
	 * Evaluate energy (store value for future use)
	 */
	public double updateEnergy() {
		if (!energyNeedsUpdate) return energy;

		energy = calcEnergy(); // Update energy

		energyNeedsUpdate = false;
		return energy;
	}

	public double[] updateGradient() {
		if (!gradientNeedsUpdate) return gradient;

		double grad[] = calcGradient();

		gradientNeedsUpdate = false;
		return grad;
	}
}
