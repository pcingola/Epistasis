package meshi.optimizers;

import meshi.energy.TotalEnergy;
import meshi.util.Terminator;

/**
 * Minimize energy according to a given set of coordinates and an energy function
 *
 **/

public abstract class Minimizer extends Optimizer {

	public final int MAX_KICKSTARTS = 1000;
	public final double tolerance;
	private double forceMagnitude;
	private int numberOfKickStrarts;
	public static final Terminator terminator = new Terminator();

	public Minimizer(TotalEnergy energy, int maxSteps, int reportEvery, double tolerance) {
		super(energy, maxSteps, reportEvery);
		this.tolerance = tolerance;
		energy.evaluate();
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
					if (testFlag) energy.test();
					throw oe;
				}
				numberOfKickStrarts++;
			}

			if (step % reportEvery == 0) System.out.println(energy().report(step));
		}
		return status(step);
	}

	private OptimizerStatus status(int step) {
		if (terminator.dead()) { return OptimizerStatus.KILLED; }
		forceMagnitude = energy.getGradMagnitude();
		if (forceMagnitude < tolerance) return OptimizerStatus.CONVERGED;
		if (step <= maxSteps) return OptimizerStatus.RUNNING;
		return OptimizerStatus.UNCONVERGED;
	}
}
