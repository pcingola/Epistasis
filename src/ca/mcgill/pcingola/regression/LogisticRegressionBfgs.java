package ca.mcgill.pcingola.regression;

import ca.mcgill.pcingola.optimizers.BFGS;

/**
 * Logistic regression
 * Model fitting by BFGS algorithm
 *
 * @author pcingola
 */
public class LogisticRegressionBfgs extends LogisticRegression {

	public LogisticRegressionBfgs(int size) {
		super(size);
		minimizer = new BFGS(this);
	}

	@Override
	public void reset() {
		super.reset();
		minimizer = new BFGS(this);
	}

}
