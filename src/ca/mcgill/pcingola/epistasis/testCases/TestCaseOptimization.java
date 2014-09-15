package ca.mcgill.pcingola.epistasis.testCases;

import java.util.Random;

import junit.framework.TestCase;
import meshi.optimizers.WolfeConditionLineSearch;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.util.Gpr;

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
		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);

		// Set initial points
		Random rand = new Random(20140915);

		int N = 100;
		for (int i = 0; i < N; i++) {
			energy.setX(0, 10.0 * rand.nextDouble() - 5.0);
			energy.setX(1, 10.0 * rand.nextDouble() - 5.0);

			double enerBef = energy.evaluate();

			// Optimize using Wolfe
			ls.findStepLength();

			double enerAfter = energy.evaluate();

			if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
			Assert.assertTrue(enerBef > enerAfter);
		}
	}

	//	public void test_03() throws Exception {
	//		TestsEnergy01 energy = new TestsEnergy01();
	//		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);
	//
	//		// Set initial value
	//		energy.setX(0, 0.8);
	//		energy.setX(1, +4);
	//
	//		double enerBef = energy.evaluate();
	//		double alpha = ls.findStepLength();
	//		double enerAfter = energy.evaluate();
	//
	//		if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
	//		Assert.assertTrue(enerBef > enerAfter);
	//	}
	//
	//	public void test_04() throws Exception {
	//		TestsEnergy01 energy = new TestsEnergy01();
	//		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);
	//
	//		// Set initial points
	//		energy.setX(0, -2);
	//		energy.setX(1, +3);
	//
	//		double enerBef = energy.evaluate();
	//		double alpha = ls.findStepLength();
	//		double enerAfter = energy.evaluate();
	//
	//		if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
	//		Assert.assertTrue(enerBef > enerAfter);
	//	}
	//
	//	public void test_05() throws Exception {
	//		TestsEnergy01 energy = new TestsEnergy01();
	//		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);
	//
	//		// Set initial points (0, 0)
	//		energy.setX(0, 0);
	//		energy.setX(1, 0);
	//
	//		double enerBef = energy.evaluate();
	//		double alpha = ls.findStepLength();
	//		double enerAfter = energy.evaluate();
	//
	//		if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
	//		Assert.assertTrue(enerBef > enerAfter);
	//	}

}
