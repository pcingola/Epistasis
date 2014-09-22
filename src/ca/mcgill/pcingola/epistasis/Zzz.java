package ca.mcgill.pcingola.epistasis;

import java.util.Random;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * Test
 *
 * @author pcingola
 */
public class Zzz {

	public static double[] realModel = { 1, -1, 0.5 };
	Random rand = new Random(20140912);
	int N = 10000;
	int size = realModel.length - 1;
	double beta[] = new double[size + 1];
	LogisticRegression lr;

	public static final boolean debug = false;

	/**
	 * Running glm in R:
	 *
	 * Call:  glm(formula = y ~ x1 + x2, family = binomial, data = dres)
	 *
	 * Coefficients:
	 * 	(Intercept)           x1           x2
	 * 	-0.6983       2.0570      -1.0052
	 *
	 * Degrees of Freedom: 199 Total (i.e. Null);  197 Residual
	 * Null Deviance:	    261.4
	 * Residual Deviance: 207.9 	AIC: 213.9
	 */
	public static void main(String[] args) {
		Timer.showStdErr("Start");

		Zzz zzz = new Zzz();
		zzz.logisticModel();
		zzz.gradient();

		Timer.showStdErr("End");
	}

	/**
	 * BGFS fitting
	 */
	public void bfgs() {
	}

	/**
	 * Gradient descent fitting
	 */
	void gradient() {
		double beta[] = new double[realModel.length];

		//		beta[0] = 2.0569899;
		//		beta[1] = -1.0051820;
		//		beta[2] = -0.6982514;

		for (int i = 0; i < realModel.length; i++)
			beta[i] = realModel[i];

		// Likelihood
		double ll = lr.logLikelihood() / Math.log(10.0);
		double llnull = lr.logLikelihoodNull() / Math.log(10.0);
		System.out.println("Log likelihood [10]: " + ll);
		System.out.println("Log likelihood Null [10]: " + llnull);

		// Learn
		lr.initModelRand();

		beta[0] = beta[1] = beta[2] = 0;
		lr.setModel(beta);
		System.out.println(lr);

		lr.setDebug(true);
		lr.learn();
		System.out.println("Model: " + lr);
	}

	public void logisticModel() {
		// Initialize model
		lr = new LogisticRegression(size);
		lr.setRand(rand);
		lr.setModel(realModel);
		lr.setThetaBest();

		// Create samples
		double in[][] = new double[N][size];
		double out[] = new double[N];
		for (int i = 0; i < N; i++) {

			// Inputs
			for (int j = 0; j < size; j++)
				in[i][j] = 2 * rand.nextDouble() - 1.0;

			// Output
			double o = lr.predict(in[i]);
			//			Gpr.debug("in: " + Gpr.toString(in[i]) + "\tout: " + o);
			out[i] = rand.nextDouble() < o ? 1.0 : 0.0;
		}
		lr.setSamples(in, out);

		// Write samples to file
		if (debug) System.out.println(lr.toStringSamples());
		String fileName = Gpr.HOME + "/logistic.txt";
		System.out.println("Model saved to: " + fileName);
		Gpr.toFile(fileName, lr.toStringSamples());

		lr.needsUpdate();
		System.out.println("Energy: " + lr.updateEnergy());
	}
}
