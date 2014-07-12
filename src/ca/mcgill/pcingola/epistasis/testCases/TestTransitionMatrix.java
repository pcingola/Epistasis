package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;

import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrixMarkov;

/**
 * Phylogenetic trees
 *
 * @author pcingola
 */
public class TestTransitionMatrix extends TestCase {

	void matrixExpTest(double q[][], double time, double expectedResult[][]) {
		TransitionMatrixMarkov Q = new TransitionMatrixMarkov(q);
		TransitionMatrixMarkov expectedExpQ = new TransitionMatrixMarkov(expectedResult);

		Q.setCheck(false);

		System.out.println("M:\n" + Q.toStringNice());
		RealMatrix expm = Q.matrix(time);
		System.out.println("exp( " + time + " * Q):\n" + (new TransitionMatrixMarkov(expm)).toStringNice());
		RealMatrix diff = expm.subtract(expectedExpQ);

		double norm = diff.getNorm();
		System.out.println("Norm: " + norm);
		Assert.assertTrue(norm < 1e-6);
	}

	public void test_01() {
		double d[][] = { //
		{ -1.0895987, 0.6663490, 0.4484602, 0.3200064 }, //
				{ 0.6663490, -1.8240385, 0.2476899, 0.9737323 }, //
				{ 0.4484602, 0.2476899, -1.6924199, 0.4895511 }, //
				{ 0.3200064, 0.9737323, 0.4895511, -1.9830556 } //
		};

		double dexpM[][] = { //
		{ 0.4844389, 0.2703457, 0.1923003, 0.2075596 }, //
				{ 0.2703457, 0.3385023, 0.1534394, 0.2488039 }, //
				{ 0.1923003, 0.1534394, 0.2607338, 0.1585244 }, //
				{ 0.2075596, 0.2488039, 0.1585244, 0.2814976 } //
		};

		matrixExpTest(d, 1.0, dexpM);
	}

	public void test_02() {
		double d[][] = { //
		{ -1.8047672, 0.6274654, 0.5668110, 0.2651656, 0.6672857 }, //
				{ 0.6274654, -1.0824480, 0.8649350, 0.8640146, 0.8659111 }, //
				{ 0.5668110, 0.8649350, -1.5335464, 0.5487009, 0.2661866 }, //
				{ 0.2651656, 0.8640146, 0.5487009, -1.8426462, 0.6265072 }, //
				{ 0.6672857, 0.8659111, 0.2661866, 0.6265072, -1.2799135 }, //
		};

		double dexpM[][] = { //
		{ 0.5274814, 0.7768704, 0.5302562, 0.4867408, 0.6174963 }, //
				{ 0.7768704, 1.3999749, 0.9022502, 0.8629877, 1.0300446 }, //
				{ 0.5302562, 0.9022502, 0.6862028, 0.5743831, 0.6575905 }, //
				{ 0.4867408, 0.8629877, 0.5743831, 0.6052957, 0.6637528 }, //
				{ 0.6174963, 1.0300446, 0.6575905, 0.6637528, 0.8988940 }, //
		};

		matrixExpTest(d, 1.23, dexpM);
	}

}
