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
	}

	@Override
	protected void init() throws OptimizerException {
		zeta = new double[logReg.getNumSamples()];
		w = new double[logReg.getNumSamples()];
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
			zeta[i] = eta[i] + (y[i] - mu[i]) / w[i]; // TODO: This might be Int or NaN!?
		}

		if (debug) {
			Gpr.debug(logReg + "\tLL: " + logReg.logLikelihood());
			Gpr.debug("\teta  (" + eta.length + "): " + Gpr.toStringHead(eta));
			Gpr.debug("\tmu   (" + mu.length + "): " + Gpr.toStringHead(mu));
			Gpr.debug("\tw    (" + w.length + "): " + Gpr.toStringHead(w));
			Gpr.debug("\tzeta (" + zeta.length + "): " + Gpr.toStringHead(zeta));
		}

		// Step II: Solve weighted least square problem
		WeightedLinearRegression wlr = new WeightedLinearRegression();
		if (!wlr.regress(zeta, logReg.getSamplesX(), w)) {
			String msg = "Cannot perform regression:" //
					+ "\n\teta  (" + eta.length + "): " + Gpr.toStringHead(eta) //
					+ "\n\tmu   (" + mu.length + "): " + Gpr.toStringHead(mu) //
					+ "\n\tw    (" + w.length + "): " + Gpr.toStringHead(w) //
					+ "\n\tzeta (" + zeta.length + "): " + Gpr.toStringHead(zeta) //
					+ "\n\tEnergy: " + energy//
					+ "\n\tIRWLS: " + this //
			;
			Gpr.debug(msg);
			return false;
		}

		// Set new coefficients
		logReg.setTheta(wlr.getCoefficients());

		return true;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append(super.toString());
		sb.append("\titerationNum: " + iterationNum);
		sb.append("\tzeta: " + Gpr.toString(zeta));
		sb.append("\tw: " + Gpr.toString(w));
		sb.append("\t" + logReg);

		return sb.toString();
	}
}
