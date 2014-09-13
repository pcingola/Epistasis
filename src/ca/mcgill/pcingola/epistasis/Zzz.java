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

		int N = 200;
		int size = 2;

		// Output
		StringBuilder sb = new StringBuilder();
		sb.append("x1\tx2\ty\n");
		Random rand = new Random(20140912);

		// Initialize model
		LogisticRegression lr = new LogisticRegression(size);
		double beta[] = new double[size + 1];
		beta[0] = 2;
		beta[1] = -1;
		beta[2] = -0.5;
		lr.setRand(rand);
		lr.setModel(beta);

		// Create samples
		double in[][] = new double[N][size];
		double out[] = new double[N];
		for (int i = 0; i < N; i++) {

			// Inputs
			for (int j = 0; j < size; j++) {
				in[i][j] = 2 * rand.nextDouble() - 1.0;
				sb.append((j > 0 ? "\t" : "") + in[i][j]);
			}

			// Output
			double o = lr.predict(in[i]);
			out[i] = rand.nextDouble() < o ? 1.0 : 0.0;
			sb.append("\t" + out[i] + "\n");
		}
		lr.setSamples(in, out);

		// Write samples to file
		System.out.println(sb);
		Gpr.toFile(Gpr.HOME + "/logistic_01.txt", sb);

		beta[0] = 2.0569899;
		beta[1] = -1.0051820;
		beta[2] = -0.6982514;
		lr.setModel(beta);

		// Likelihood
		double ll = lr.logLikelihood() / Math.log(10.0);
		double llnull = lr.logLikelihoodNull() / Math.log(10.0);
		System.out.println("Log likelihood [10]: " + ll);
		System.out.println("Log likelihood Null [10]: " + llnull);

		// Learn
		lr.initModelRand();
		lr.setDebug(true);
		double betaModel[] = lr.learn();
		System.out.println("Model: " + lr);

		Timer.showStdErr("End");
	}
}
