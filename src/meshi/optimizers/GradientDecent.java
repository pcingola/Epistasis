package meshi.optimizers;

import meshi.optimizers.exceptions.OptimizerException;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Calculate gradient and move in that direction
 *
 */

public class GradientDecent extends Minimizer {

	double mu = 1.0; // Learning rate 

	public GradientDecent(Energy energy) {
		super(energy);
	}

	@Override
	protected void init() throws OptimizerException {
		energy().evaluate();
	}

	@Override
	protected void kickStart() throws OptimizerException {
	}

	@Override
	protected boolean minimizationStep() throws OptimizerException {
		energy.setThetaBest();
		energy.addThetaBestGradient(-mu);
		energy.evaluate();

		if (debug) Gpr.debug("Energy: " + energy.getEnergy() //
				+ "\tTheta: " + Gpr.toString(energy.getTheta()) //
				+ "\tGradient: " + Gpr.toString(energy.getTheta()) //
		);
		return true;
	}
}
