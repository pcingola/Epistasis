package ca.mcgill.pcingola.epistasis.testCases;

import java.util.Arrays;
import java.util.Random;

import junit.framework.TestCase;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.junit.Assert;

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
public class TestCaseLaplaceIntegral extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	void compareMaxDiff(double dexp[], double dreal[], double maxAbsDiff) {
		for (int i = 0; i < dexp.length; i++) {
			double err = Math.abs(dexp[i] - dreal[i]);
			if (err > maxAbsDiff) throw new RuntimeException("Error: " + err + "\n\tIndex " + i + "\tExpecting: " + dexp[i] + "\n\tObtained: " + dreal[i]);
		}
	}

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
		lr.setSamplesAddIntercept(in, out);

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
	public LogisticRegression modelFitTest(Random rand, double beta[], int N, String createFile, String loadFile, double betaFit[], double maxDifference, String minType) {
		int size = beta.length - 1;

		// Initialize model
		LogisticRegression lr = new LogisticRegression(size);
		lr.setDebug(debug);
		if (rand != null) lr.setRand(rand);
		lr.setModel(beta);

		// Create samples
		if (loadFile == null) createSamples(lr, N, size, createFile, rand);
		else readModel(lr, loadFile);

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

		compareMaxDiff(betaModel, betaFit, maxDifference);

		return lr;
	}

	/**
	 * Read a model
	 */
	public void readModel(LogisticRegression lr, String fileName) {
		String lines[] = Gpr.readFile(fileName).split("\n");

		// Model size and sampels
		int N = lines.length - 1; // Number of samples. First line is "title" (not sample)
		int size = lines[0].split("\t").length - 2; // Number of elements in linear model. First column is 'output', other columns are input (one input is 'bias')
		if (debug) Gpr.debug("File '" + fileName + "', samples: " + N + ", model size: " + size);

		double in[][] = new double[N][size];
		double out[] = new double[N];

		// Load model: Read and parse file
		int i = 0;
		for (String l : lines) {
			// Skip title
			if (i == 0) {
				i++;
				continue;
			}

			String f[] = l.split("\t");

			// Output
			out[i - 1] = Gpr.parseDoubleSafe(f[0]);

			// Inputs (first two columns are 'output' and 'intercept')
			for (int j = 2; j < f.length; j++)
				in[i - 1][j - 2] = Gpr.parseDoubleSafe(f[j]);

			i++;
		}

		// Show matrix
		if (debug) Gpr.debug("In:\n" + Gpr.toString(in));

		// Set samples
		lr.setSamplesAddIntercept(in, out);
	}

	public void test_01_laplace() {
		Gpr.debug("Test");
		Random rand = new Random(20140912);
		int N = 300;

		double beta[] = { 2, -1 }; // Real model
		double betaFit[] = { 2.2634, -1.037 }; // Expected fitted model

		//---
		// Calculate integral using Laplace's approximation
		//---
		LogisticRegression lr = modelFitTest(rand, beta, N, null, null, betaFit, 0.01, "irwls");
		if (debug) Gpr.debug("Beta (model fit  ): " + Gpr.toString(lr.getTheta()) + "\tLogLik : " + lr.logLikelihood());

		int n = beta.length;
		double H[][] = new double[n][n];
		double sx[][] = lr.getSamplesX(); // Samples input values
		lr.getSamplesY();
		double p[] = lr.getOut(); // Output of logistic regresion is probability

		// Calculate the Hessian matrix
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				double ss = 0;
				for (int s = 0; s < N; s++) {
					ss += sx[s][i] * sx[s][j] * p[s] * (1 - p[s]);
				}

				H[i][j] = ss;
			}

		if (debug) Gpr.debug("H:\n" + Gpr.toString(H));

		// Calculate Hessian's determinant
		Array2DRowRealMatrix Hr = new Array2DRowRealMatrix(H);
		double detH = (new LUDecomposition(Hr)).getDeterminant();
		if (debug) Gpr.debug("det(H): " + detH);

		// Use Lapplace formula
		double likelihood = Math.exp(lr.logLikelihood());
		double intLaplace = likelihood * 2.0 * Math.PI * Math.sqrt(1.0 / detH);
		if (debug) Gpr.debug("Integral (Laplace): " + intLaplace);

		//---
		// Calculate integral using "brute force" approach
		//---
		double step = 0.05;
		double xmin = -5;
		double xmax = 5;

		double llmax = Double.NEGATIVE_INFINITY;
		double betaBfMax[] = null;
		double sum = 0;
		for (double x0 = xmin; x0 < xmax; x0 += step)
			for (double x1 = xmin; x1 < xmax; x1 += step) {
				// Set parameters
				lr.setTheta(0, x0);
				lr.setTheta(1, x1);

				// Calculate likelihood and add
				double ll = lr.logLikelihood();
				sum += Math.exp(ll);

				// Is this the maximum log likelihood so far?
				if (ll > llmax) {
					llmax = ll;
					betaBfMax = Arrays.copyOf(lr.getTheta(), beta.length);
				}
			}

		// Show and check maximum likelihood
		if (debug) Gpr.debug("Beta (brute force): " + Gpr.toString(betaBfMax) + "\tLogLik : " + llmax);
		compareMaxDiff(betaFit, betaBfMax, step);

		// Integral
		double integral = sum * (step * step);
		if (debug) Gpr.debug("Integral (brute force): " + integral);

		// Check that integrals have simmilar values
		double diff = Math.abs((integral - intLaplace) / integral);
		Assert.assertTrue(diff < 0.01);

	}
}
