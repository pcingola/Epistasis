package meshi.optimizers;

import meshi.optimizers.exceptions.OptimizerException;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Minimize energy according to a given set of coordinates and an energy function
 *
 **/

public abstract class Minimizer extends Optimizer {

	public final int MAX_KICKSTARTS = 100;

	int numberOfKickStrarts;

	public Minimizer(Energy energy) {
		super(energy);
	}

	/**
	 * Initialize algorithm
	 */
	protected abstract void init() throws OptimizerException;

	/**
	 * Kick-start the algorithm after minimization step failed
	 */
	protected abstract void kickStart() throws OptimizerException;

	/**
	 * Minimize
	 * @return true if OK, false if 'kickstart' is required
	 * @throws OptimizerException
	 */
	protected abstract boolean minimizationStep() throws OptimizerException;

	@Override
	public OptimizerStatus run() {
		try {
			return run(true);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	public OptimizerStatus run(boolean testFlag) throws OptimizerException {
		init();
		numberOfKickStrarts = 0;
		int step;
		for (step = 1; optimizerTerminator.status(step) == OptimizerStatus.RUNNING; step++) {

			energy.evaluate(); // Update energy

			boolean minimizationStepOK = minimizationStep();

			if (debug) Gpr.debug("Iteration: " + step + "\t" + this);

			if (!minimizationStepOK) {
				if (numberOfKickStrarts >= MAX_KICKSTARTS) {
					optimizerTerminator.kill("\n\nThe simulation was restarted for " + MAX_KICKSTARTS + " times " //
							+ "which is more than allowed.\n" //
							+ "So many restarts are indicative of an ill-shaped energy function or " //
							+ "an energy differentiation\n" //
							);
				} else {
					try {
						kickStart();
						if (verbose) System.err.println("kickstart # " + numberOfKickStrarts + " done");
					} catch (OptimizerException oe) {
						throw oe;
					}
					numberOfKickStrarts++;
				}
			}
		}

		if (debug) Gpr.debug(this);
		return optimizerTerminator.status(step);
	}

	@Override
	public String toString() {
		return this.getClass().getSimpleName() + " : " + energy;
	}

}
