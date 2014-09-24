package ca.mcgill.pcingola.epistasis.testCases;

import java.util.Random;

import junit.framework.TestCase;
import meshi.optimizers.BFGS;
import meshi.optimizers.GradientDecent;
import meshi.optimizers.Minimizer;
import meshi.optimizers.SteepestDecent;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseLogisticRegression extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = true || debug;

	/**
	 * Create a model, create data, fit model and check that parameters are correctly fitted
	 * @param rand : Random number generator
	 * @param beta : Model to create
	 * @param N : Number of data samples
	 * @param createFile : If non-null, create a tab-separated file with sample data
	 */
	public void modelFitTest(Random rand, double beta[], int N, String createFile, double betaFit[], double maxDifference, String minType) {
		int size = beta.length - 1;

		// Output file titles
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < size; i++)
			sb.append((i > 0 ? "\t" : "") + "x" + i);
		sb.append("\ty\n");

		// Initialize model
		LogisticRegression lr = new LogisticRegression(size);
		lr.setDebug(debug);
		lr.setRand(rand);
		lr.setModel(beta);

		// Minimizer
		Minimizer minimizer = null;
		switch (minType) {
		case "bfgs":
			minimizer = new BFGS(lr);
			break;

		case "steepest":
			minimizer = new SteepestDecent(lr);
			break;

		case "grad":
			minimizer = new GradientDecent(lr);
			break;

		default:
			throw new RuntimeException("UNknown type " + minType);
		}
		lr.setMinnimizer(minimizer);

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
		if (createFile != null) {
			System.out.println(sb);
			Gpr.toFile(createFile, sb);
		}

		// Learn

		lr.initModelRand();
		double betaModel[] = lr.learn();
		if (verbose) System.out.println("Model " + minType + ": " + lr);

		for (int i = 0; i < beta.length; i++) {
			double err = Math.abs(betaModel[i] - betaFit[i]);
			if (err > maxDifference) throw new RuntimeException("Error: " + err + "\n\tExpecting: " + betaFit[i] + "\n\tObtained: " + betaModel[i]);
		}
	}

	public void test_01() {
		Random rand = new Random(20140912);
		int N = 200;

		double beta[] = { 2, -1, -0.5 }; // Real model
		double betaFit[] = { 2.055550258008242, -1.0041789502014213, -0.6979724967536511 }; // Expected fitted model

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "grad");
	}

	public void test_01_bfgs() {
		Random rand = new Random(20140912);
		int N = 200;

		double beta[] = { 2, -1, -0.5 }; // Real model
		double betaFit[] = { 2.057, -1.005, -0.698 }; // Expected fitted model

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "bfgs");
	}

	public void test_01_steepest() {
		Random rand = new Random(20140912);
		int N = 200;

		double beta[] = { 2, -1, -0.5 }; // Real model
		double betaFit[] = { 2.057, -1.005, -0.698 }; // Expected fitted model

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "steepest");
	}

	public void test_02() {
		Random rand = new Random(20140912);
		int N = 10000;

		double beta[] = { 2, -1, -0.5 }; // Real model
		double betaFit[] = { 2.02, -0.9848, -0.5247 }; // Expected fitted model

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "grad");
	}

	public void test_02_bfgs() {
		Random rand = new Random(20140912);
		int N = 100000;

		double beta[] = { 2, -1, -0.5 }; // Real model
		double betaFit[] = { 2.015, -1.015, -0.506 }; // Expected fitted model

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "bfgs");
	}

	public void test_02_steepest() {
		Random rand = new Random(20140912);
		int N = 100000;

		double beta[] = { 2, -1, -0.5 }; // Real model
		double betaFit[] = { 2.015, -1.015, -0.506 }; // Expected fitted model

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "steepest");
	}

	public void test_03() {
		Random rand = new Random(20140912);
		int N = 10000;
		double beta[] = { 1.7, -0.1, -1, 0.8, 1.3, -0.5 };
		double betaFit[] = { 1.6997, -0.0427, -1.0012, 0.8036, 1.2781, -0.5192 };

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "grad");
	}

}
