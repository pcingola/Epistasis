package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;
import meshi.optimizers.SteepestDecent;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseOptimization extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = true;

	public void test_01() {
		TestsEnergy01 energy = new TestsEnergy01();
		SteepestDecent optimizer = new SteepestDecent(energy);
		optimizer.run();

	}

}
