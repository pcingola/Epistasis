package meshi.optimizers;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A generic class for all energy functions to be minimized.
 */
public abstract class Energy {

	protected double[] x;
	protected double[] xBest;
	protected double[] gradient;
	protected double energy;
	protected boolean energyNeedsUpdate = true;
	protected boolean gradientNeedsUpdate = true;

	public Energy(int dim) {
		x = new double[dim];
		xBest = new double[dim];
		gradient = new double[dim];
		energy = Double.NaN;
	}

	public void addX(int idx, double addToXi) {
		x[idx] += addToXi;
		energyNeedsUpdate = gradientNeedsUpdate = true;
	}

	public void addXBestGradient(double alpha) {
		for (int i = 0; i < x.length; i++)
			x[i] = xBest[i] + alpha * gradient[i];

		energyNeedsUpdate = gradientNeedsUpdate = true;
	}

	public void addXGradient(double alpha) {
		for (int i = 0; i < x.length; i++)
			x[i] += alpha * gradient[i];

		energyNeedsUpdate = gradientNeedsUpdate = true;
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
	 * Copy current x[] to an array
	 */
	public void copyX(double copyX[]) {
		System.arraycopy(x, 0, copyX, 0, x.length);
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

	public double[] getX() {
		return x;
	}

	public void setX(double newX[]) {
		System.arraycopy(newX, 0, x, 0, x.length);
		energyNeedsUpdate = gradientNeedsUpdate = true;
	}

	public void setX(int idx, double xi) {
		x[idx] = xi;
		energyNeedsUpdate = gradientNeedsUpdate = true;
	}

	public void setXBest() {
		for (int i = 0; i < x.length; i++)
			x[i] = xBest[i];

		energyNeedsUpdate = gradientNeedsUpdate = true;
	}

	public int size() {
		return x.length;
	}

	@Override
	public String toString() {
		return "Energy: " + (energyNeedsUpdate ? "[Needs update]" : "") //
				+ energy //
				//
				+ "\tx: " + Gpr.toString(x) //
				//
				+ "\tgradient " //
				+ (gradientNeedsUpdate ? "[Needs update]" : "") //
				+ ": " + Gpr.toString(gradient);
	}

	/**
	 * Evaluate energy (store value for future use)
	 */
	public double updateEnergy() {
		if (!energyNeedsUpdate) return energy;

		double energyOld = energy;
		energy = calcEnergy(); // Update energy

		// Best energy so far? Keep a copy
		if (energy < energyOld) System.arraycopy(x, 0, xBest, 0, x.length);

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
