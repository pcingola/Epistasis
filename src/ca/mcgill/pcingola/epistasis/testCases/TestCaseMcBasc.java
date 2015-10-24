package ca.mcgill.pcingola.epistasis.testCases;

import org.junit.Assert;
import org.junit.Test;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.coEvolutionMetrics.AaSimilarityMatrix;
import ca.mcgill.pcingola.epistasis.coEvolutionMetrics.McBasc;
import junit.framework.TestCase;

/**
 * Test cases for McBASC algorithm
 *
 * @author pcingola
 */
public class TestCaseMcBasc extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	@Test
	public void test01() {
		Gpr.debug("Test");

		// Load matrix
		String matrix = "data/McLachlan_matrix.txt";
		AaSimilarityMatrix simMatrix = new AaSimilarityMatrix(matrix);

		// Create an 'alignment'
		String coli = "AGH";
		String colj = "AGH";

		// Invoke McBASB algorithm
		McBasc mcBasc = new McBasc(simMatrix, coli, colj);
		mcBasc.setDebug(debug);
		double score = mcBasc.score();
		if (verbose) System.out.println("Score: " + mcBasc.score());

		// Check results
		Assert.assertEquals(1.4444444444444444, mcBasc.getMeanS_i(), 1e-6);
		Assert.assertEquals(1.4444444444444444, mcBasc.getMeanS_j(), 1e-6);

		Assert.assertEquals(2.697735676039774, mcBasc.getSigmaS_i(), 1e-6);
		Assert.assertEquals(2.697735676039774, mcBasc.getSigmaS_j(), 1e-6);

		Assert.assertEquals(score, 0.30986714363098145, 1e-6);
	}

}
