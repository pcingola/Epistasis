package ca.mcgill.pcingola.epistasis.testCases;

import java.util.Random;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.optimizers.BFGS;
import ca.mcgill.pcingola.optimizers.GradientDecent;
import ca.mcgill.pcingola.optimizers.IRWLS;
import ca.mcgill.pcingola.optimizers.Minimizer;
import ca.mcgill.pcingola.optimizers.SteepestDecent;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseLogisticRegression extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	/**
	 * Create samples and set them to logReg model (optionally create a file)
	 */
	void createSamples(LogisticRegression lr, int N, int size, String createFile, Random rand) {
		// Output file titles
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < size; i++)
			sb.append((i > 0 ? "\t" : "") + "x" + i);
		sb.append("\ty\n");

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

	}

	/**
	 * Create a model, create data, fit model and check that parameters are correctly fitted
	 * @param rand : Random number generator
	 * @param beta : Model to create
	 * @param N : Number of data samples
	 * @param createFile : If non-null, create a tab-separated file with sample data
	 */
	public LogisticRegression modelFitTest(Random rand, double beta[], int N, String createFile, double betaFit[], double maxDifference, String minType) {
		int size = beta.length - 1;
		// Initialize model
		LogisticRegression lr = new LogisticRegression(size);
		lr.setDebug(debug);
		if (rand != null) lr.setRand(rand);
		lr.setModel(beta);

		// Create samples
		createSamples(lr, N, size, createFile, rand);

		// Minimizer
		Minimizer minimizer = null;
		switch (minType) {
		case "bfgs":
			minimizer = new BFGS(lr);
			break;

		case "irwls":
			minimizer = new IRWLS(lr);
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

		// Learn
		double betaModel[] = lr.learn();
		if (verbose) System.out.println("Model " + minType + ": " + lr);

		for (int i = 0; i < beta.length; i++) {
			double err = Math.abs(betaModel[i] - betaFit[i]);
			if (err > maxDifference) throw new RuntimeException("Error: " + err + "\n\tExpecting: " + betaFit[i] + "\n\tObtained: " + betaModel[i]);
		}

		return lr;
	}

	//
	//	public void test_01() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 200;
	//
	//		double beta[] = { 2, -1, -0.5 }; // Real model
	//		double betaFit[] = { 2.055550258008242, -1.0041789502014213, -0.6979724967536511 }; // Expected fitted model
	//
	//		modelFitTest(rand, beta, N, null, betaFit, 0.01, "grad");
	//	}
	//
	//	public void test_01_bfgs() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 200;
	//
	//		double beta[] = { 2, -1, -0.5 }; // Real model
	//		double betaFit[] = { 2.057, -1.005, -0.698 }; // Expected fitted model
	//
	//		modelFitTest(rand, beta, N, null, betaFit, 0.01, "bfgs");
	//	}
	//
	//	public void test_01_steepest() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 200;
	//
	//		double beta[] = { 2, -1, -0.5 }; // Real model
	//		double betaFit[] = { 2.057, -1.005, -0.698 }; // Expected fitted model
	//
	//		modelFitTest(rand, beta, N, null, betaFit, 0.01, "steepest");
	//	}
	//
	//	public void test_02() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 10000;
	//
	//		double beta[] = { 2, -1, -0.5 }; // Real model
	//		double betaFit[] = { 2.02, -0.9848, -0.5247 }; // Expected fitted model
	//
	//		modelFitTest(rand, beta, N, null, betaFit, 0.01, "grad");
	//	}
	//
	//	public void test_02_bfgs() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 100000;
	//
	//		double beta[] = { 2, -1, -0.5 }; // Real model
	//		double betaFit[] = { 2.015, -1.015, -0.506 }; // Expected fitted model
	//
	//		modelFitTest(rand, beta, N, null, betaFit, 0.01, "bfgs");
	//	}
	//
	//	public void test_02_steepest() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 100000;
	//
	//		double beta[] = { 2, -1, -0.5 }; // Real model
	//		double betaFit[] = { 2.015, -1.015, -0.506 }; // Expected fitted model
	//
	//		modelFitTest(rand, beta, N, null, betaFit, 0.01, "steepest");
	//	}
	//
	//	public void test_03() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 10000;
	//		double beta[] = { 1.7, -0.1, -1, 0.8, 1.3, -0.5 };
	//		double betaFit[] = { 1.6997, -0.0427, -1.0012, 0.8036, 1.2781, -0.5192 };
	//
	//		modelFitTest(rand, beta, N, null, betaFit, 0.01, "grad");
	//	}
	//
	//	/**
	//	 * Reset model and learn same data
	//	 */
	//	public void test_04() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 10000;
	//		double beta[] = { 1.7, -0.1, -1, 0.8, 1.3, -0.5 };
	//		double betaFit[] = { 1.6997, -0.0427, -1.0012, 0.8036, 1.2781, -0.5192 };
	//
	//		LogisticRegression lr = modelFitTest(rand, beta, N, null, betaFit, 0.01, "bfgs");
	//
	//		// Make sure we can learn after a reset
	//		lr.reset();
	//		rand = new Random(20140912);
	//		double betaScond[] = { 1.7, -0.1, -1, 0.8, 1.3, -0.5 };
	//		double betaFitSecond[] = { 1.6997, -0.0427, -1.0012, 0.8036, 1.2781, -0.5192 };
	//		lr = modelFitTest(rand, betaScond, N, null, betaFitSecond, 0.01, "bfgs");
	//	}
	//
	//	/**
	//	 * Reset model and learn different data
	//	 */
	//	public void test_05() {
	//		Gpr.debug("Test");
	//		Random rand = new Random(20140912);
	//		int N = 10000;
	//		double beta[] = { 1.7, -0.1, -1, 0.8, 1.3, -0.5 };
	//		double betaFit[] = { 1.6997, -0.0427, -1.0012, 0.8036, 1.2781, -0.5192 };
	//
	//		LogisticRegression lr = modelFitTest(rand, beta, N, null, betaFit, 0.01, "bfgs");
	//
	//		// Make sure we can learn after a reset
	//		lr.reset();
	//		rand = new Random(20140912);
	//		double betaScond[] = { -0.7, 0.8, 2, -0.8, 1.9, 0.5 };
	//		double betaFitSecond[] = { -0.705, 0.834, 1.980, -0.772, 1.889, 0.531 };
	//		lr = modelFitTest(rand, betaScond, N, null, betaFitSecond, 0.01, "bfgs");
	//	}
	//
	//	public void test_06() {
	//		throw new RuntimeException("Create test using T2D data. Compare to R's calculations");
	//	}

	public void test_01_irwls() {
		Gpr.debug("Test");
		Random rand = new Random(20140912);
		int N = 50;

		double beta[] = { 2, -1, -0.5 }; // Real model
		double betaFit[] = { 2.057, -1.005, -0.698 }; // Expected fitted model

		modelFitTest(rand, beta, N, null, betaFit, 0.01, "irwls");
	}

}
