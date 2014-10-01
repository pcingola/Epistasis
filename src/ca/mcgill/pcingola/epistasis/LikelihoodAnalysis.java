package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.StreamSupport;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.regression.LogisticRegression;
import ca.mcgill.pcingola.regression.LogisticRegressionBfgs;

/**
 * Logistic regression log-likelihood analysis of VCF + phenotype data
 *
 * @author pcingola
 */
public class LikelihoodAnalysis {

	public static boolean writeToFile = false;

	public static final int PHENO_ROW_NUMBER = 0; // Covariate number zero is phenotype

	String phenoCovariatesFileName = Gpr.HOME + "/t2d1/coEvolution/coEvolution.pheno.covariates.txt";
	String vcfFileName = Gpr.HOME + "/t2d1/vcf/eff/hm.chr1.gt.vcf";
	//	String vcfFileName = Gpr.HOME + "/t2d1/vcf/eff/t2d_13K.test_02.vcf";

	boolean debug = false;
	int numSamples, numCovs;
	int count = 0;
	int covariatesToNormalize[] = { 11, 12 };
	double covariates[][];
	double pheno[];
	double llMax = Double.NEGATIVE_INFINITY, llMin = Double.POSITIVE_INFINITY;
	String logLikInfoField; // If not null, an INFO field is added
	String sampleIds[];
	LogisticRegression lr;
	HashMap<Long, LogisticRegression> modelAltByThread = new HashMap<Long, LogisticRegression>();
	HashMap<Long, LogisticRegression> modelNullByThread = new HashMap<Long, LogisticRegression>();

	public static void main(String[] args) {
		Timer.showStdErr("Start");

		boolean debug = false;

		LikelihoodAnalysis zzz = new LikelihoodAnalysis(args);

		if (debug) {
			zzz.setDebug(debug);
			zzz.setLogLikInfoField("LL");
		}

		zzz.run(debug);

		Timer.showStdErr("End");
	}

	public LikelihoodAnalysis(String args[]) {
		if (args.length > 1) {
			phenoCovariatesFileName = args[0];
			vcfFileName = args[1];
		}
	}

	/**
	 * Create models
	 */
	void createModels() {
		Timer.showStdErr("Creating models for thread: " + Thread.currentThread());
		long threadId = Thread.currentThread().getId();

		//---
		// Create alternative model
		//---
		LogisticRegression lrAlt = new LogisticRegressionBfgs(numCovs);

		// Copy all covariates, except first one (phenotype row)
		double xAlt[][] = new double[numSamples][numCovs];
		for (int i = 0; i < numSamples; i++)
			for (int j = 0; j < numCovs; j++)
				xAlt[i][j] = covariates[i][j];

		lrAlt.setSamples(xAlt, pheno);
		modelAltByThread.put(threadId, lrAlt);

		//---
		// Create null model
		//---
		LogisticRegression lrNull = new LogisticRegressionBfgs(numCovs - 1); // No genotypes

		// Copy all covariates, except first one (phenotype row)
		double xNull[][] = new double[numSamples][numCovs - 1];
		for (int i = 0; i < numSamples; i++)
			for (int j = 0; j < numCovs - 1; j++)
				xNull[i][j] = covariates[i][j + 1];

		lrNull.setSamples(xNull, pheno);
		modelNullByThread.put(threadId, lrNull);

		// Samples to skip
		boolean skip[] = new boolean[numSamples];
		Arrays.fill(skip, false);
		lrAlt.setSkip(skip);
		lrNull.setSkip(skip);
	}

	/**
	 * Load phenotypes and covariates
	 */
	void loadPhenoAndCovariates(String phenoCovariates) {
		Timer.showStdErr("Reading PCA data form '" + phenoCovariates + "'");

		// Read "phenotypes + covariates" file
		String lines[] = Gpr.readFile(phenoCovariates).split("\n");

		// Allocate PCA matrix
		numCovs = lines.length - 1; // Row 1 is the title, other rows are covariates
		numSamples = lines[0].split("\t").length - 1; // All, but first column are samples
		covariates = new double[numSamples][numCovs];

		// Parse
		int covariateNum = -1;
		for (String line : lines) {
			line = line.trim();
			String fields[] = line.split("\t");

			// Skip title
			if (covariateNum < 0) {
				sampleIds = new String[fields.length - 1];
				for (int i = 1; i < fields.length; i++)
					sampleIds[i - 1] = fields[i];
			} else {
				// Parse covariate values
				for (int i = 1; i < fields.length; i++)
					covariates[i - 1][covariateNum] = Gpr.parseDoubleSafe(fields[i]);
			}

			covariateNum++;
		}

		//---
		// Convert phenotype to {0,1}
		//---
		pheno = new double[numSamples];
		for (int i = 0; i < numSamples; i++)
			pheno[i] = covariates[i][PHENO_ROW_NUMBER] - 1.0;

		//---
		// Normalize covariates: sex, age
		//---
		for (int cov : covariatesToNormalize)
			normalizeCovariates(cov);
	}

	/**
	 * Fit models and calculate log likelihood test
	 */
	void logLikelihood(VcfEntry ve) {
		// Get models for this thread
		long threadId = Thread.currentThread().getId();
		LogisticRegression lrAlt = modelAltByThread.get(threadId);
		LogisticRegression lrNull = modelNullByThread.get(threadId);

		// Need to create models?
		if (lrAlt == null) {
			createModels();
			lrAlt = modelAltByThread.get(threadId);
			lrNull = modelNullByThread.get(threadId);
		}

		// Reset models
		lrAlt.reset();
		lrNull.reset();

		// Get genotypes
		byte gt[] = ve.getGenotypesScores();

		// Copy genotypes to first row (alt model) and set 'skip' field
		double xAlt[][] = lrAlt.getSamplesX();
		boolean skip[] = lrAlt.getSkip();
		for (int vcfSampleNum = 0; vcfSampleNum < numSamples; vcfSampleNum++) {
			xAlt[vcfSampleNum][PHENO_ROW_NUMBER] = gt[vcfSampleNum];
			skip[vcfSampleNum] = (gt[vcfSampleNum] < 0) || (pheno[vcfSampleNum] < 0);
		}

		// Calculate logistic models
		lrAlt.setDebug(debug);
		lrAlt.learn();

		lrNull.setDebug(debug);
		lrNull.learn();

		// Calculate likelihood ratio
		double ll = -2.0 * (lrAlt.logLikelihood() - lrNull.logLikelihood());

		if (logLikInfoField != null) ve.addInfo(logLikInfoField, "" + ll);

		// Stats
		if (Double.isFinite(ll)) {
			boolean show = (llMax < ll) || (llMin > ll);
			show = true;

			llMin = Math.min(llMin, ll);
			llMax = Math.max(llMax, ll);

			if (show) {
				System.out.println(ve.toStr() //
						+ "\tLL_ratio: " + ll //
						+ "\tLL range: " + llMin + " / " + llMax //
						+ "\n\tModel Alt  : " + lrAlt //
						+ "\n\tModel Null : " + lrNull //
				);
			} else Timer.show(count + "\t" + ve.toStr());

		} else {
			throw new RuntimeException("Likelihood ratio is infinite!\n" + ve);
		}

		// TODO: Calculate and check p-value (Chi-square test)

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

	public void run() {
		run(false);
	}

	/**
	 * Run analysis and collect some results (for test cases)
	 */
	public List<VcfEntry> run(boolean createList) {
		//---
		// Load phenotype and covariates
		//---
		loadPhenoAndCovariates(phenoCovariatesFileName);

		//---
		// Read VCF file and perform logistic regression
		// TODO: Matrix format
		//			- Use "SnpSift allelmat" format it's much faster
		//			- We need to load the whole matrix into memory (filter out singletons / MAF < 0.1% ?)
		//			- Perform a quick check for overlapping variants between two positions before using logistic regression
		//---
		VcfFileIterator vcf = new VcfFileIterator(vcfFileName);

		// TODO
		Gpr.debug("WRITE TEST CASE TO COMPARE TO R's RESULTS");

		//---
		// Check that sample names and sample order matches
		//---
		vcf.readHeader();
		List<String> sampleNames = vcf.getVcfHeader().getSampleNames();
		int snum = 0;
		for (String s : sampleNames) {
			if (!s.equals(sampleIds[snum])) { throw new RuntimeException("Sample names do not match:" //
					+ "\n\tSample [" + snum + "] in VCF file        :  '" + s + "'" //
					+ "\n\tSample [" + snum + "] in phenotypes file :  '" + sampleIds[snum] + "'" //
			); }
			snum++;
		}

		//---
		// Calculate for each entry in VCF file (use parallel stream
		//---
		ArrayList<VcfEntry> list = new ArrayList<VcfEntry>();
		if (createList) {
			// Process and populate list of VCF entries (single thread, used for debugging and test cases)
			for (VcfEntry ve : vcf) {
				logLikelihood(ve);
				list.add(ve);
			}
		} else StreamSupport.stream(vcf.spliterator(), true).forEach(ve -> logLikelihood(ve)); // Process (do not populate list)

		return list;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setLogLikInfoField(String logLikInfoField) {
		this.logLikInfoField = logLikInfoField;
	}
}
