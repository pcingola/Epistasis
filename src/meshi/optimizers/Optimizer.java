package meshi.optimizers;

import meshi.energy.TotalEnergy;
import meshi.util.Terminator;

/**
 * Minimize energy according to a given set of coordinates and an energy function
 *
 */
public abstract class Optimizer {
	public enum OptimizerStatus {
		RUNNING, CONVERGED, UNCONVERGED, KILLED, DONE;
	}

	public final TotalEnergy energy;
	public final int maxSteps;
	public final int reportEvery;
	public static final Terminator optimizerTerminator = new Terminator();

	public Optimizer(TotalEnergy energy, int maxSteps, int reportEvery) {
		this.maxSteps = maxSteps;
		this.energy = energy;
		this.reportEvery = reportEvery;
		optimizerTerminator.reset();
	}

	public TotalEnergy energy() {
		return energy;
	}

	public abstract OptimizerStatus run() throws OptimizerException;
}
