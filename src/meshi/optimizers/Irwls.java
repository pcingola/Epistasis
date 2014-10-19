package meshi.optimizers;

import meshi.optimizers.exceptions.OptimizerException;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * IRWLS optimization algorithm (specifically implemented for Logistic regression)
 * References:
 * 
 * 	i) "Introduction to Generalized Linear Models" by Heather Turner
 *		http://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf
 *
 *	ii) "Exercises (Part 4), Introduction to R UCLA/CCPR", John Fox, February 2005
 *		http://socserv.socsci.mcmaster.ca/jfox/Courses/UCLA/Exercises-4.pdf
 *
 *															Pablo Cingolani 2014
 **/

public class Irwls extends Minimizer {

	int iterationNum = 0;
	LogisticRegression logReg;
	double nu[]; // Variance parameter
	double zeta[]; // Output and derivate
	double w[]; // Weights for re-weighted least squares

	public Irwls(LogisticRegression logReg) {
		super(logReg);
		this.logReg = logReg;
		nu = new double[logReg.getNumSamples()];
		zeta = new double[logReg.getNumSamples()];
		w = new double[logReg.getNumSamples()];
	}

	@Override
	protected void init() throws OptimizerException {
		Gpr.debug("Init");
	}

	/**
	 * Starting the BFGS minimization by a few steepest descent steps, followed by inverse Hessian initialization
	 */
	@Override
	protected void kickStart() throws OptimizerException {
		if (debug) Gpr.debug("A kick start has occurred in iteration:" + iterationNum + "\n");
	}

	/**
	 * IRWLS algorithm
	 */
	@Override
	protected boolean minimizationStep() throws OptimizerException {
		logReg.evaluate();

		logReg.getH();
		logReg.getOut();

		return true;
	}
}
