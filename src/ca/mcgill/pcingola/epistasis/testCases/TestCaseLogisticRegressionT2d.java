package ca.mcgill.pcingola.epistasis.testCases;

import java.util.List;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.LikelihoodAnalysis;

/**
 * Test cases for logistic regression using phenotypes and VCF data
 *
 * @author pcingola
 */
public class TestCaseLogisticRegressionT2d extends TestCase {

	public static boolean debug = true;
	public static boolean verbose = false || debug;

	/**
	 * Reset model and learn different data
	 */
	public void test_01() {
		Gpr.debug("Test");

		String args[] = { "test/pheno.covariates.T2D_13K.txt", "test/t2d_13K.test_00.vcf" };
		LikelihoodAnalysis la = new LikelihoodAnalysis(args);
		la.setDebug(debug);

		String llInfo = "LL";
		la.setLogLikInfoField(llInfo);

		List<VcfEntry> list = la.run(true);

		// Check result (only on line)
		System.out.println(list.get(0).getInfo(llInfo));
	}
}
