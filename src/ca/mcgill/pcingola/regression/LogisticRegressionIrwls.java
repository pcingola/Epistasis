package ca.mcgill.pcingola.regression;

import ca.mcgill.pcingola.optimizers.IRWLS;

/**
 * Logistic regression
 * Model fitting by IRWLS algorithm
 *
 * @author pcingola
 */
public class LogisticRegressionIrwls extends LogisticRegression {

	public LogisticRegressionIrwls(int size) {
		super(size);
		minimizer = new IRWLS(this);
	}

	@Override
	public void reset() {
		super.reset();
		minimizer = new IRWLS(this);
	}

}
