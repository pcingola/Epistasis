package ca.mcgill.pcingola.epistasis;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import meshi.optimizers.BFGS;
import meshi.optimizers.GradientDecent;
import meshi.optimizers.Minimizer;
import meshi.optimizers.SteepestDecent;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.regression.LogisticRegression;
import ca.mcgill.pcingola.regression.LogisticRegressionBfgs;

/**
 * Test
 *
 * @author pcingola
 */
public class Zzz {

	public static double[] realModel = { 2, -1, -0.5 };

	public static final boolean debug = false;
	public static String type = "grad"; // "steepest"; "bgfs";

	String t2dPhenoCovariates = Gpr.HOME + "/t2d1/coEvolution/coEvolution.pheno.covariates.txt";
	String t2dVcf = Gpr.HOME + "/t2d1/vcf/eff/hm.chr1.gt.vcf";
	// String t2dVcf = Gpr.HOME + "/t2d1/vcf/eff/hm.test.vcf";
	HashMap<String, Integer> sampleId2pos;
	Random rand = new Random(20140912);
	int N = 10000;
	int size = realModel.length - 1;
	double beta[] = new double[size + 1];
	double covariates[][];
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
		zzz.logisticT2d();

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
		double ll = lr.logLikelihood();
		System.out.println("Log likelihood: " + ll);

		// Learn
		lr.initModelRand();

		beta[0] = beta[1] = beta[2] = 0;
		lr.setModel(beta);
		System.out.println(lr);

		lr.setDebug(true);
		lr.learn();
	}

	/**
	 * Load phenotypes and covariates
	 */
	void loadPhenoAndCovariates(String phenoCovariates) {
		Timer.showStdErr("Reading PCA data form '" + phenoCovariates + "'");

		// Read "phenotypes + covariates" file
		String lines[] = Gpr.readFile(phenoCovariates).split("\n");
		sampleId2pos = new HashMap<>();

		// Allocate PCA matrix
		int numCov = lines.length - 1; // Row 1 is the title, other rows are covariates
		int numSamples = lines[0].split("\t").length - 1; // All, but first column are samples
		covariates = new double[numSamples][numCov];

		// Parse
		int covariateNum = -1;
		String sampleIds[];
		for (String line : lines) {
			line = line.trim();
			String fields[] = line.split("\t");

			// Skip title
			if (covariateNum < 0) {
				sampleIds = fields;

				// Add sampleIDs to map
				for (int i = 1; i < fields.length; i++)
					sampleId2pos.put(sampleIds[i], i - 1);

			} else {
				// Parse covariate values
				for (int i = 1; i < fields.length; i++)
					covariates[i - 1][covariateNum] = Gpr.parseDoubleSafe(fields[i]);
			}

			covariateNum++;
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

	/***
	 * Logistic regression using T2D-26K data
	 *
	 * Test dataset:
	 * 		- VCF          : t2d1/vcf/eff/hm.chr1.gt.vcf
	 * 		- pheno + PCAs : t2d1/coEvolution/coEvolution.pheno.covariates.txt
	 * 		- R script     : workspace/Epistasis/scripts/coEvolution/coEvolution.r
	 */
	void logisticT2d() {
		boolean writeToFile = false;

		//---
		// Load phenotype and covariates
		//---
		loadPhenoAndCovariates(t2dPhenoCovariates);

		int numSamples = covariates.length;
		int numCovs = covariates[0].length;

		// TODO: Convert phenotype to {0,1}
		int phenoRowNum = 0; // Covariate number zero is phenotype
		double pheno[] = new double[numSamples];
		for (int i = 0; i < numSamples; i++)
			pheno[i] = covariates[i][phenoRowNum] - 1.0;

		// Normalize covariates: sex, age
		normalizeCovariates(11);
		normalizeCovariates(12);

		//---
		// Create logistic regression models
		//---

		// Initialize logistic regression
		// Note: Even thought the first row currently has the
		//       phenotype, we'll overwrite the row using 'allele'
		//       information before performing logistic regression

		LogisticRegression lrAlt = new LogisticRegressionBfgs(numCovs);
		lrAlt.setSamples(covariates, pheno);

		LogisticRegression lrNull = new LogisticRegressionBfgs(numCovs - 1); // No genotypes

		// Copy all covariates, except first one (phenotype row)
		double xNull[][] = new double[numSamples][numCovs - 1];
		for (int i = 0; i < numSamples; i++)
			for (int j = 0; j < numCovs - 1; j++)
				xNull[i][j] = covariates[i][j + 1];

		lrNull.setSamples(xNull, pheno);

		// Samples to skip
		boolean skip[] = new boolean[numSamples];
		Arrays.fill(skip, false);
		lrAlt.setSkip(skip);
		lrNull.setSkip(skip);

		//---
		// Read VCF file and perform logistic regression
		// TODO: Matrix format
		//			- Use "SnpSift allelmat" format it's much faster
		//			- We need to load the whole matrix into memory (filter out singletons / MAF < 0.1% ?)
		//			- Perform a quick check for overlapping variants between two positions before using logistic regression
		//---
		VcfFileIterator vcf = new VcfFileIterator(t2dVcf);

		double llMax = Double.NEGATIVE_INFINITY, llMin = Double.POSITIVE_INFINITY;
		// TODO: Parallelize using
		//       StreamSupport.stream(vcf.spliterator(), true);

		int count = 1;
		for (VcfEntry ve : vcf) {

			// Reset model
			lrAlt.reset();
			lrNull.reset();

			// Get genotypes
			byte gt[] = ve.getGenotypesScores();

			// Copy genotypes to first row and set 'skip' field
			double xAlt[][] = lrAlt.getSamplesX();
			for (int i = 0; i < numSamples; i++) {
				xAlt[i][phenoRowNum] = gt[i];
				skip[i] = (gt[i] < 0) || (pheno[i] < 0);
			}

			// Calculate logistic models
			lrAlt.setDebug(debug);
			lrAlt.learn();

			lrNull.setDebug(debug);
			lrNull.learn();

			// Calculate likelihood ratio
			double ll = 2.0 * (lrAlt.logLikelihood() - lrNull.logLikelihood());

			//
			if (Double.isFinite(ll)) {
				if ((llMax < ll) || (llMin > ll)) {
					System.out.println(ve.toStr() //
							+ "\tLL_ratio: " + ll //
							+ "\tLL range: " + llMin + " / " + llMax //
							+ "\n\tModel Alt  : " + lrAlt //
							+ "\n\tModel Null : " + lrNull //
					);
				} else Timer.show(count + "\t" + ve.toStr());

				llMin = Math.min(llMin, ll);
				llMax = Math.max(llMax, ll);
			} else {
				throw new RuntimeException("Likelihood ratio is infinite!\n" + ve);
			}

			// TODO: Calculate and check p-value
			//         - Chi-square
			//         - Wald test

			// Save as TXT table
			if (writeToFile) {
				String fileName = Gpr.HOME + "/lr_test.alt.txt";
				Gpr.debug("Writing 'alt' table to :" + fileName);
				Gpr.toFile(fileName, lrAlt.toStringSamples());

				fileName = Gpr.HOME + "/lr_test.null.txt";
				Gpr.debug("Writing 'null' table to :" + fileName);
				Gpr.toFile(fileName, lrNull.toStringSamples());
			}

			count++;
		}

		// TODO
		Gpr.debug("WRITE TEST CASE TO COMPARE TO R's RESULTS");
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
		double ll = lr.logLikelihood();
		System.out.println("Log likelihood: " + ll);

		Timer.showStdErr("End");
	}

	/**
	 * Normalize (mean 0, variance 1) a covariates' row
	 */
	void normalizeCovariates(int covNum) {
		// Calculate mean
		double sum = 0;
		for (int i = 0; i < covariates.length; i++)
			sum += covariates[i][covNum];
		double avg = sum / covariates.length;

		// Calculate stddev
		sum = 0;
		for (int i = 0; i < covariates.length; i++) {
			double d = (covariates[i][covNum] - avg);
			sum += d * d;
		}
		double stddev = Math.sqrt(sum / (covariates.length - 1));

		Timer.showStdErr("Covariate " + covNum + ", mean: " + avg + ", stddev: " + stddev);

		// Normalize
		for (int i = 0; i < covariates.length; i++)
			covariates[i][covNum] = (covariates[i][covNum] - avg) / stddev;

	}
}
