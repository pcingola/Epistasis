package ca.mcgill.pcingola.epistasis.testCases;

import java.util.Random;

import junit.framework.TestCase;
import meshi.optimizers.SteepestDecent;
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

	public static boolean debug = false;
	public static boolean verbose = false;

	public void test_01() {
		TestsEnergy01 energy = new TestsEnergy01();
		SteepestDecent optimizer = new SteepestDecent(energy);
		optimizer.setDebug(debug);
		optimizer.run();

		Assert.assertEquals(energy.getTheta()[0], 2.0, EPSILON);
		Assert.assertEquals(energy.getTheta()[1], 5.0, EPSILON);
	}

	public void test_02() throws Exception {
		TestsEnergy01 energy = new TestsEnergy01();
		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);

		// Set initial points
		Random rand = new Random(20140915);

		int N = 10000;
		for (int i = 0; i < N; i++) {
			energy.setTheta(0, 10.0 * rand.nextDouble() - 5.0);
			energy.setTheta(1, 10.0 * rand.nextDouble() - 5.0);
			if (debug) System.out.println("\n\nx: " + Gpr.toString(energy.getTheta()));

			double enerBef = energy.evaluate();

			// Optimize using Wolfe
			ls.reset();
			ls.setDebug(debug);
			ls.findStepLength();

			double enerAfter = energy.evaluate();

			if (verbose) Gpr.debug("Energy before: " + enerBef + "\tafter: " + energy);
			Assert.assertTrue(enerBef > enerAfter);
		}
	}

	public void test_03() throws Exception {
		TestsEnergy01 energy = new TestsEnergy01();
		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);
		ls.setDebug(debug);

		// Set initial value
		energy.setTheta(0, 0.8);
		energy.setTheta(1, +4);

		double enerBef = energy.evaluate();
		double alpha = ls.findStepLength();
		double enerAfter = energy.evaluate();

		if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
		Assert.assertTrue(enerBef > enerAfter);
	}

	public void test_04() throws Exception {
		TestsEnergy01 energy = new TestsEnergy01();
		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);

		// Set initial points
		energy.setTheta(0, -2);
		energy.setTheta(1, +3);

		double enerBef = energy.evaluate();
		double alpha = ls.findStepLength();
		double enerAfter = energy.evaluate();

		if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
		Assert.assertTrue(enerBef > enerAfter);
	}

	public void test_05() throws Exception {
		TestsEnergy01 energy = new TestsEnergy01();
		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);

		// Set initial points (0, 0)
		energy.setTheta(0, 0);
		energy.setTheta(1, 0);

		double enerBef = energy.evaluate();
		double alpha = ls.findStepLength();
		double enerAfter = energy.evaluate();

		if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
		Assert.assertTrue(enerBef > enerAfter);
	}

	public void test_06() throws Exception {
		TestsEnergy01 energy = new TestsEnergy01();
		WolfeConditionLineSearch ls = new WolfeConditionLineSearch(energy);

		// Set initial points (0, 0)
		energy.setTheta(0, 2.784);
		energy.setTheta(1, -2.108);

		double enerBef = energy.evaluate();
		double alpha = ls.findStepLength();
		double enerAfter = energy.evaluate();

		if (debug) Gpr.debug("Energy before: " + enerBef + "\tafter: " + enerAfter + "\n\t" + energy);
		Assert.assertTrue(enerBef > enerAfter);
	}

}
