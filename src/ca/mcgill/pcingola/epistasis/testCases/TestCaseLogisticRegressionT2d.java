package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.LikelihoodAnalysis;

/**
 * Test cases for logistic regression using phenotypes and VCF data
 *
 * @author pcingola
 */
public class TestCaseLogisticRegressionT2d extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	/**
	 * Reset model and learn different data
	 */
	public void test_01() {
		Gpr.debug("Test");

		String args[] = { "test/pheno.covariates.T2D_13K.txt", "test/t2d_13K.test_01.vcf" };
		LikelihoodAnalysis la = new LikelihoodAnalysis(args);
		la.run();
	}

}
