package meshi.energy;

/**
 * A generic class for all energy functions to be minimized.
 */
public abstract class Energy {

	double[][] coordinates;
	double energy;

	/**
	 * Calculate energy
	 */
	public abstract double calcEnergy();

	public double[][] coordinates() {
		return coordinates;
	}

	/**
	 * Evaluate energy (store value for future use)
	 */
	public double evaluate() {
		energy = calcEnergy();
		return energy;
	}

	/**
	 * Returns latest energy value
	 **/
	public double getEnergy() {
		return energy;
	}
}
