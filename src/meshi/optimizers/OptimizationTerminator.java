package meshi.optimizers;

import meshi.optimizers.Optimizer.OptimizerStatus;

/**
 * Control when the optimization step is finished
 *
 * @author pcingola
 */
public class OptimizationTerminator {

	public static final int DEFAULT_MAX_STEPS = 1000000;
	public static final double DEFAULT_GRADIENT_MAX_ABS_THRESHOLD = 1E-6;

	boolean dead;
	int maxSteps;
	double gradientMaxAbsThreshold;
	double energyOld = Double.MAX_VALUE;
	String message;
	Energy energy;

	public OptimizationTerminator(Energy energy) {
		this(energy, DEFAULT_MAX_STEPS, DEFAULT_GRADIENT_MAX_ABS_THRESHOLD);
	}

	public OptimizationTerminator(Energy energy, int maxSteps, double gradientMaxAbsThreshold) {
		this.energy = energy;
		this.maxSteps = maxSteps;
		this.gradientMaxAbsThreshold = gradientMaxAbsThreshold;
		message = "";
		dead = false;
	}

	/**
	 * Finds the maximal component (in magnitude) of the gradient vector in coordinates ( coordinates[][1] ).
	 */
	double getGradMagnitude() {
		double grad[] = energy.updateGradient();

		double maxAbs = 0;
		for (int i = 0; i < grad.length; i++)
			maxAbs = Math.max(Math.abs(grad[i]), maxAbs);

		return maxAbs;
	}

	public boolean isDead() {
		return dead;
	}

	public void kill(String message) {
		this.message = message;
		dead = true;
	}

	public String message() {
		return message;
	}

	public void setGradientMaxAbsThreshold(double gradientMaxAbsThreshold) {
		this.gradientMaxAbsThreshold = gradientMaxAbsThreshold;
	}

	public void setMaxSteps(int maxSteps) {
		this.maxSteps = maxSteps;
	}

	public OptimizerStatus status(int step) {
		// Has this been killed?
		if (isDead()) { return OptimizerStatus.KILLED; }

		// Is energy improving?
		double energyNew = energy.getEnergy();
		if (energyNew >= energyOld) return OptimizerStatus.CONVERGED;
		energyOld = energyNew;

		// Is the gradient 'strong' enough
		double gradMaxAbs = getGradMagnitude();
		if (gradMaxAbs < gradientMaxAbsThreshold) return OptimizerStatus.CONVERGED;

		// Are we done with the number of steps?
		if (step <= maxSteps) return OptimizerStatus.RUNNING;

		return OptimizerStatus.UNCONVERGED;
	}

}
