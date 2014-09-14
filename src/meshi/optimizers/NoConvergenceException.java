package meshi.optimizers;
import java.util.*;

public class NoConvergenceException extends Exception {	
	public NoConvergenceException(int nSteps) {
		super("Minimization did not converge after "+nSteps+" steps");		
	}
}

