package ca.mcgill.pcingola.regression;

import meshi.optimizers.GradientDecent;

/**
 * Logistic regression
 * Model fitting by BFGS algorithm
 *
 * @author pcingola
 */
public class LogisticRegressionGrad extends LogisticRegression {

	public LogisticRegressionGrad(int size) {
		super(size);
		minimizer = new GradientDecent(this);
	}

	@Override
	public void reset() {
		super.reset();
		minimizer = new GradientDecent(this);
	}

}
