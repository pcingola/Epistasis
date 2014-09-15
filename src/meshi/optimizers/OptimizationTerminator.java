package meshi.optimizers;

import meshi.optimizers.Optimizer.OptimizerStatus;

public class OptimizationTerminator {

	public static final int DEFAULT_MAX_STEPS = 1000000;
	public static final double DEFAULT_GRADIENT_MAX_ABS_THRESHOLD = 1E-6;

	boolean dead;
	int maxSteps;
	double gradientMaxAbsThreshold;
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
		double grad[] = energy.getGradient();

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

	public OptimizerStatus status(int step) {
		if (isDead()) { return OptimizerStatus.KILLED; }

		double gradMaxAbs = getGradMagnitude();
		if (gradMaxAbs < gradientMaxAbsThreshold) return OptimizerStatus.CONVERGED;

		if (step <= maxSteps) return OptimizerStatus.RUNNING;

		return OptimizerStatus.UNCONVERGED;
	}

}
