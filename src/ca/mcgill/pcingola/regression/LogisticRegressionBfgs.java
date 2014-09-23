package ca.mcgill.pcingola.regression;

import meshi.optimizers.BFGS;

/**
 * Logistic regression
 * Model fitting by BFGS algorithm
 *
 * @author pcingola
 */
public class LogisticRegressionBfgs extends LogisticRegression {

	public LogisticRegressionBfgs(int size) {
		super(size);
		minnimizer = new BFGS(this);
	}

}
