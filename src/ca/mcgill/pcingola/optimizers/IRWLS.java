package ca.mcgill.pcingola.optimizers;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.optimizers.exceptions.OptimizerException;
import ca.mcgill.pcingola.regression.LogisticRegression;
import ca.mcgill.pcingola.regression.WeightedLinearRegression;

/**
 * IRWLS optimization algorithm (Iterated Re-Weighted Least Squares)
 * Specifically implemented for Logistic Regression (binary output)
 *
 * References:
 *
 * 	i) "Introduction to Generalized Linear Models" by Heather Turner
 *		http://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf
 *
 *	ii) "Exercises (Part 4), Introduction to R UCLA/CCPR", John Fox, February 2005
 *		http://socserv.socsci.mcmaster.ca/jfox/Courses/UCLA/Exercises-4.pdf
 *
 * Note: We use "Generalized Linear models" nomenclature
 *															Pablo Cingolani 2014
 **/

public class IRWLS extends Minimizer {

	int iterationNum = 0;
	LogisticRegression logReg;
	double zeta[]; // Output and derivate
	double w[]; // Weights for re-weighted least squares

	public IRWLS(LogisticRegression logReg) {
		super(logReg);
		this.logReg = logReg;
		zeta = new double[logReg.getNumSamples()];
		w = new double[logReg.getNumSamples()];
	}

	@Override
	protected void init() throws OptimizerException {
		// Initialize beta to 0
		double zero[] = new double[energy.getDim()];
		energy.setTheta(zero);
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
		// Step I: Evaluate logistic regression and
		//         calculate intermediate variables nu, zeta, w
		logReg.evaluate();
		double eta[] = logReg.getH();
		double mu[] = logReg.getOut();
		double y[] = logReg.getSamplesY();

		int n = logReg.getNumSamples();
		for (int i = 0; i < n; i++) {
			w[i] = mu[i] * (1.0 - mu[i]);
			zeta[i] = eta[i] + (y[i] - mu[i]) / w[i];
		}

		if (debug) {
			Gpr.debug(logReg + "\tLL: " + logReg.logLikelihood());
			Gpr.debug("\teta  : " + Gpr.toString(eta));
			Gpr.debug("\tmu   : " + Gpr.toString(mu));
			Gpr.debug("\tw    : " + Gpr.toString(w));
			Gpr.debug("\tzeta : " + Gpr.toString(zeta));
		}

		// Step II: Solve weighted least square problem
		WeightedLinearRegression wlr = new WeightedLinearRegression();
		if (!wlr.regress(zeta, logReg.getSamplesX(), w)) return false;

		// Set new coefficients
		logReg.setTheta(wlr.getCoefficients());

		return true;
	}
}
