package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.GwasEpistasis;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseGwas extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	public void test_01_Gwas_Map() {
		Gpr.debug("Test");

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String genesLikeFile = "test/NM_001438.txt";
		String vcfFile = ""; // It doesn't matter, it is not used

		GwasEpistasis gwasEpistasis = new GwasEpistasis(configFile, genome, genesLikeFile, vcfFile);
		gwasEpistasis.setDebug(debug);
		gwasEpistasis.initialize();
		gwasEpistasis.readGenesLogLikelihood();

		// Check that all mappings are OK
		Assert.assertEquals(0, gwasEpistasis.getCountErr());
		Assert.assertEquals(334, gwasEpistasis.getCountOk());
	}

	public void test_02_Gwas_Map() {
		Gpr.debug("Test");

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String genesLikeFile = "test/NM_021969.txt";
		String vcfFile = ""; // It doesn't matter, it is not used

		GwasEpistasis gwasEpistasis = new GwasEpistasis(configFile, genome, genesLikeFile, vcfFile);
		gwasEpistasis.setDebug(debug);
		gwasEpistasis.initialize();
		gwasEpistasis.readGenesLogLikelihood();

		// Check that all mappings are OK
		Assert.assertEquals(0, gwasEpistasis.getCountErr());
		Assert.assertEquals(458, gwasEpistasis.getCountOk());
	}
}
