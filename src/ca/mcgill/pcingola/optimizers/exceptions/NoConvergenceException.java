package ca.mcgill.pcingola.optimizers.exceptions;

public class NoConvergenceException extends Exception {

	public NoConvergenceException(int nSteps) {
		super("Minimization did not converge after " + nSteps + " steps");
	}
}
