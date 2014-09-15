package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;
import meshi.optimizers.WolfConditionLineSearch;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseOptimization extends TestCase {

	public static final double EPSILON = 0.0000001;

	public static boolean debug = true;
	public static boolean verbose = true;

	//	public void test_01() {
	//		TestsEnergy01 energy = new TestsEnergy01();
	//		SteepestDecent optimizer = new SteepestDecent(energy);
	//		optimizer.setDebug(debug);
	//		optimizer.run();
	//
	//		Assert.assertEquals(energy.getX()[0], 2.0, EPSILON);
	//		Assert.assertEquals(energy.getX()[1], 5.0, EPSILON);
	//	}

	public void test_02() throws Exception {
		TestsEnergy01 energy = new TestsEnergy01();
		WolfConditionLineSearch ls = new WolfConditionLineSearch(energy);

		energy.evaluate();
		ls.findStepLength();
	}

}
