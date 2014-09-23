package meshi.optimizers;

import meshi.optimizers.exceptions.OptimizerException;

/**
 * Minimize energy according to a given set of coordinates and an energy function
 */
public abstract class Optimizer {

	public enum OptimizerStatus {
		RUNNING, CONVERGED, UNCONVERGED, KILLED, DONE;
	}

	protected boolean debug = false;
	protected boolean verbose = false || debug;
	protected Energy energy;
	protected OptimizationTerminator optimizerTerminator;

	public Optimizer(Energy energy) {
		this.energy = energy;
		optimizerTerminator = new OptimizationTerminator(energy);
	}

	public Energy energy() {
		return energy;
	}

	public abstract OptimizerStatus run() throws OptimizerException;

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	@Override
	public String toString() {
		return this.getClass().getSimpleName() + " : " + energy;
	}

}
