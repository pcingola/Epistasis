package meshi.optimizers;

import meshi.energy.Energy;
import meshi.optimizers.exceptions.OptimizerException;

/**
 * Minimize energy according to a given set of coordinates and an energy function
 *
 **/

public abstract class Minimizer extends Optimizer {

	public final int MAX_KICKSTARTS = 1000;
	public final double tolerance;
	private double forceMagnitude;
	private int numberOfKickStrarts;
	public static final OptimizationTerminator terminator = new OptimizationTerminator();

	public Minimizer(Energy energy, int maxSteps, int reportEvery, double tolerance) {
		super(energy, maxSteps, reportEvery);
		this.tolerance = tolerance;
		energy.evaluate();
	}

	/**
	 * Finds the maximal component (in magnitude) of the gradient vector in coordinates ( coordinates[][1] ).
	 */
	double getGradMagnitude() {
		throw new RuntimeException("Gradient magnitude: Unimplemented!");
	}

	protected abstract void init() throws OptimizerException;

	protected abstract void kickStart() throws OptimizerException;

	protected abstract boolean minimizationStep() throws OptimizerException;

	@Override
	public OptimizerStatus run() throws OptimizerException {
		return run(true);
	}

	public OptimizerStatus run(boolean testFlag) throws OptimizerException {
		init();
		numberOfKickStrarts = 0;
		int step;
		for (step = 1; status(step) == OptimizerStatus.RUNNING; step++) {
			boolean minimizationStepOK = minimizationStep();

			if (!minimizationStepOK) {
				if (numberOfKickStrarts >= MAX_KICKSTARTS) throw new OptimizerException("\n\nThe simulation was restarted for " + MAX_KICKSTARTS + " times " + "which is more than allowed.\n" + "So many restarts are indicative of an ill-shaped energy function or " + "an energy differentiation\n");
				try {
					kickStart();
					System.out.println("kickstart # " + numberOfKickStrarts + " done");
				} catch (OptimizerException oe) {
					throw oe;
				}
				numberOfKickStrarts++;
			}
		}

		return status(step);
	}

	private OptimizerStatus status(int step) {
		if (terminator.dead()) { return OptimizerStatus.KILLED; }
		forceMagnitude = getGradMagnitude();
		if (forceMagnitude < tolerance) return OptimizerStatus.CONVERGED;
		if (step <= maxSteps) return OptimizerStatus.RUNNING;
		return OptimizerStatus.UNCONVERGED;
	}
}
