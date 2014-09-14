package meshi.optimizers;

import meshi.energy.Energy;
import meshi.optimizers.exceptions.OptimizerException;

/**
 * Minimize energy according to a given set of coordinates and an energy function
 */
public abstract class Optimizer {

	public enum OptimizerStatus {
		RUNNING, CONVERGED, UNCONVERGED, KILLED, DONE;
	}

	public final Energy energy;
	public final int maxSteps;
	public final int reportEvery;
	public static final OptimizationTerminator optimizerTerminator = new OptimizationTerminator();

	public Optimizer(Energy energy, int maxSteps, int reportEvery) {
		this.maxSteps = maxSteps;
		this.energy = energy;
		this.reportEvery = reportEvery;
		optimizerTerminator.reset();
	}

	public Energy energy() {
		return energy;
	}

	public abstract OptimizerStatus run() throws OptimizerException;
}
