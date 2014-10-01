package ca.mcgill.pcingola.epistasis.testCases;

import java.util.List;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.LikelihoodAnalysis;

/**
 * Test cases for logistic regression using phenotypes and VCF data
 *
 * @author pcingola
 */
public class TestCaseLogisticRegressionT2d extends TestCase {

	public static final String LL_INFO_FIELD = "LL";

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	/**
	 * Model: Genotype is 0/0 for all (so it should match the null model)
	 */
	public void test_01() {
		Timer.showStdErr("Test");

		String args[] = { "test/pheno.covariates.T2D_13K.txt", "test/t2d_13K.test_01.vcf" };
		LikelihoodAnalysis la = new LikelihoodAnalysis(args);
		la.setDebug(debug);
		la.setWriteToFile(debug);
		la.setLogLikInfoField(LL_INFO_FIELD);

		List<VcfEntry> list = la.run(true);

		// Check result (only on line)
		System.out.println(list.get(0).getInfo(LL_INFO_FIELD));
		Timer.showStdErr("done.");
		throw new RuntimeException("Missing check condition!");
	}

	//	/**
	//	 * Model: Genotype is phenotype (so it should match perfect)
	//	 */
	//	public void test_02() {
	//		Gpr.debug("Test");
	//
	//		String args[] = { "test/pheno.covariates.T2D_13K.txt", "test/t2d_13K.test_02.vcf" };
	//		LikelihoodAnalysis la = new LikelihoodAnalysis(args);
	//		la.setWriteToFile(debug);
	//		la.setDebug(debug);
	//		la.setLogLikInfoField(LL_INFO_FIELD);
	//
	//		List<VcfEntry> list = la.run(true);
	//
	//		// Check result (only on line)
	//		System.out.println(list.get(0).getInfo(LL_INFO_FIELD));
	//		throw new RuntimeException("Missing check condition!");
	//	}

}
