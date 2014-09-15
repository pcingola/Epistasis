package ca.mcgill.pcingola.epistasis.testCases;

import meshi.optimizers.Energy;

/**
 * Implements an "Energy" function:
 *
 * 		(x - 2)^2 + (y - 5)^2
 *
 * @author pcingola
 */

public class TestsEnergy01 extends Energy {

	public TestsEnergy01() {
		super(2);
	}

	@Override
	protected double calcEnergy() {
		double d1 = (x[0] - 2.0);
		double d2 = (x[1] - 5.0);
		return (d1 * d1) + (d2 * d2);
	}

	@Override
	protected double[] calcGradient() {
		gradient[0] = 2.0 * (x[0] - 2.0);
		gradient[1] = 2.0 * (x[1] - 5.0);
		return gradient;
	}

}
