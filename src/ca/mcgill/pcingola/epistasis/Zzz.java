package ca.mcgill.pcingola.epistasis;

import java.util.HashMap;
import java.util.Random;

import meshi.optimizers.BFGS;
import meshi.optimizers.GradientDecent;
import meshi.optimizers.Minimizer;
import meshi.optimizers.SteepestDecent;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.regression.LogisticRegression;

/**
 * Test
 *
 * @author pcingola
 */
public class Zzz {

	public static double[] realModel = { 2, -1, -0.5 };

	public static final boolean debug = false;
	public static String type = "grad"; // "steepest"; "bgfs";

	String eigen26k = "data/26k/pruned_v3_based_variants_26k_pca.txt.gz";
	HashMap<String, Integer> sampleId2pos;
	Random rand = new Random(20140912);
	int N = 10000;
	int size = realModel.length - 1;
	double beta[] = new double[size + 1];
	double pca[][];
	LogisticRegression lr;

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
		// zzz.logisticTest();
		zzz.logistic26k();

		Timer.showStdErr("End");
	}

	/**
	 * Learn: Fit model
	 */
	public void learn() {
		double beta[] = new double[realModel.length];

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
	}

	/***
	 * Logistic regression using T2D-26K data
	 */
	void logistic26k() {
		Timer.showStdErr("Reading PCA data form '" + eigen26k + "'");
		// Read file
		String lines[] = Gpr.readFile(eigen26k).split("\n");
		sampleId2pos = new HashMap<>();

		// Allocate PCA matrix
		int numPca = lines[0].split("\t").length - 2; // All, but sampleId and population are PCAs
		pca = new double[lines.length][numPca];

		// Parse
		int sampleNum = 0;
		for (String line : lines) {
			line = line.trim();
			String fields[] = line.split("\t");

			// Parse sample ID & population data
			String sampleId = fields[0];
			String population = fields[fields.length - 1];

			// Add to map
			sampleId2pos.put(sampleId, sampleNum);

			// Parse PCA values
			if (debug) System.out.print(sampleId + "\t" + population);
			for (int i = 1; i < fields.length - 1; i++) {
				pca[sampleNum][i - 1] = Gpr.parseDoubleSafe(fields[i]);
				if (debug) System.out.print("\t" + pca[sampleNum][i - 1]);
			}
			if (debug) System.out.println("");

			sampleNum++;
		}
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

	void logisticTest() {
		// Create model
		logisticModel();

		// Select minimizer type and learn
		Minimizer minimizer = null;
		switch (type) {
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
			throw new RuntimeException("UNknown type " + type);
		}
		lr.setMinnimizer(minimizer);

		learn();

		// Show model after fitting
		System.out.println("Model: " + lr);
		double ll = lr.logLikelihood() / Math.log(10.0);
		double llnull = lr.logLikelihoodNull() / Math.log(10.0);
		System.out.println("Log likelihood [10]: " + ll);
		System.out.println("Log likelihood Null [10]: " + llnull);

		Timer.showStdErr("End");
	}
}
