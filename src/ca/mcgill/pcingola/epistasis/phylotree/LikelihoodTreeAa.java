package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;

import org.apache.commons.math3.linear.RealMatrix;

import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * A phylogenetic tree that can be used to calculate likelihoods
 *
 * Ref: Felsestein's algorithm
 * 		"Computational Molecular Evolution" Z. Yang, Section 4.2.2, page 102
 * 		"Biological Sequence Analysis", Durbin et. al., page 199
 *
 * @author pcingola
 */
public class LikelihoodTreeAa extends PhylogeneticTree {

	public static final double GAP_PROB = 1.0 / (GprSeq.AMINO_ACIDS.length);
	double p[];

	/**
	 * Create root node
	 */
	public LikelihoodTreeAa() {
		super();
	}

	public LikelihoodTreeAa(PhylogeneticTree parent, String phyloStr) {
		super(parent, phyloStr);
	}

	public LikelihoodTreeAa(String name) {
		super(name);
	}

	public LikelihoodTreeAa(String name, PhylogeneticTree left, double distanceLeft, PhylogeneticTree right, double distanceRight) {
		super(name, left, distanceLeft, right, distanceRight);
	}

	public double[] getP() {
		return p;
	}

	public double getP(int seqCode) {
		return p[seqCode];
	}

	/**
	 * Calculate likelihood
	 */
	public double likelihood(TransitionMatrix tmatrix, double pi[]) {
		double likelihood = 0.0;
		resetNode(pi.length);

		for (int aaCode = 0; aaCode < p.length; aaCode++) {
			p[aaCode] = likelihood(tmatrix, aaCode);
			likelihood += p[aaCode] * pi[aaCode];
		}

		return likelihood;
	}

	/**
	 * Calculate likelihood for this seqCode
	 */
	protected double likelihood(TransitionMatrix tmatrix, int aaCode) {
		// Already calculated?
		if (!Double.isNaN(p[aaCode])) return p[aaCode];

		//---
		// Leaf node?
		//---
		if (isLeaf()) {
			if (aaCode < 0) return GAP_PROB; // Uniform probability for GAPs

			// Probability is 1 for that sequence, 0 for others
			if (sequenceCode < 0) p[aaCode] = GAP_PROB; // Uniform probability for GAPs
			else if (sequenceCode == aaCode) p[aaCode] = 1.0;
			else p[aaCode] = 0.0;
			return p[aaCode];
		}

		//---
		// Non-leaf node
		//---

		// Likelihood from the left sub-tree
		RealMatrix P = tmatrix.matrix(distanceLeft);
		double pleft = 0;
		if (left != null) {
			for (int sc = 0; sc < p.length; sc++) {
				double lleft = ((LikelihoodTreeAa) left).likelihood(tmatrix, sc);
				double pij = P.getEntry(aaCode, sc);
				pleft += lleft * pij;
			}
		} else pleft = 1.0; // No node

		// Likelihood from the right sub-tree
		P = tmatrix.matrix(distanceRight);
		double pright = 0;
		if (right != null) {
			for (int sc = 0; sc < p.length; sc++) {
				double lright = ((LikelihoodTreeAa) right).likelihood(tmatrix, sc);
				double pij = P.getEntry(aaCode, sc);
				pright += lright * pij;
			}
		} else pright = 1.0; // No node

		// Set code
		p[aaCode] = pleft * pright;

		return p[aaCode];
	}

	@Override
	protected PhylogeneticTree newNode(PhylogeneticTree parent, String phyloStr) {
		return new LikelihoodTreeAa(parent, phyloStr);
	}

	/**
	 * Pre-calculate matrix exponentials
	 */
	public void precalculateExpm(TransitionMatrix tmatrix) {
		if (isLeaf()) return;

		Timer.showStdErr("Dim(Q): " + tmatrix.getRowDimension() + "\tExp(" + distanceLeft + ")");
		tmatrix.matrix(distanceLeft);
		((LikelihoodTreeAa) left).precalculateExpm(tmatrix);

		Timer.showStdErr("Dim(Q): " + tmatrix.getRowDimension() + "\tExp(" + distanceRight + ")");
		tmatrix.matrix(distanceRight);
		((LikelihoodTreeAa) right).precalculateExpm(tmatrix);
	}

	@Override
	protected void resetNode(int size) {
		if (p == null) p = new double[size];
		Arrays.fill(p, Double.NaN);

		if (left != null) left.resetNode(size);
		if (right != null) right.resetNode(size);
	}

	public String toStringP() {
		StringBuilder sb = new StringBuilder();

		sb.append(name + " '" + getSequence() + "' :\t[ ");
		sb.append(p[0]);
		for (int i = 1; i < p.length; i++)
			sb.append(", " + p[i]);
		sb.append(" ]");

		return sb.toString();
	}

}
