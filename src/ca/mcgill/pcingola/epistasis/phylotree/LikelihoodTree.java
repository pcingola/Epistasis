package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;

import org.apache.commons.math3.linear.RealMatrix;

import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * A phylogenetic tree that can be used to calculate likelihoods
 *
 * Ref: Felsestein's algorithm
 * 		"Computational Molecular Evolution" Z. Yang, Section 4.2.2, page 102
 * 		"Biological Sequence Analysis", Durbin et. al., page 199
 *
 * @author pcingola
 */
public class LikelihoodTree extends PhylogeneticTree {

	public static final double GEP_PROB = 1.0 / (GprSeq.AMINO_ACIDS.length);
	double p[];

	/**
	 * Create root node
	 */
	public LikelihoodTree() {
		super();
		resetNode();
	}

	public LikelihoodTree(PhylogeneticTree parent, String phyloStr) {
		super(parent, phyloStr);
		resetNode();
	}

	public LikelihoodTree(String name) {
		super(name);
		resetNode();
	}

	public LikelihoodTree(String name, PhylogeneticTree left, double distanceLeft, PhylogeneticTree right, double distanceRight) {
		super(name, left, distanceLeft, right, distanceRight);
		resetNode();
	}

	public double[] getP() {
		return p;
	}

	public double getP(int seqCode) {
		return p[seqCode];
	}

	/**
	 * Calculate likelihood
	 * @param tmatrix : Transition matrix
	 * @return
	 */
	public double likelihood(TransitionMatrix tmatrix, double pi[]) {
		double likelihood = 0.0;
		reset();

		for (int scode = 0; scode < p.length; scode++) {
			p[scode] = likelihood(tmatrix, scode);
			likelihood += p[scode] * pi[scode];
		}

		return likelihood;
	}

	/**
	 * Calculate likelihood for this seqCode
	 * @param tmatrix
	 * @param seqCode
	 * @return
	 */
	protected double likelihood(TransitionMatrix tmatrix, int seqCode) {
		// Already calculated?
		if (!Double.isNaN(p[seqCode])) return p[seqCode];

		//---
		// Leaf node?
		//---
		if (isLeaf()) {
			if (seqCode < 0) return GEP_PROB; // Uniform probability for GAPs

			// Probability is 1 for that sequence, 0 for others
			if (sequenceCode < 0) p[seqCode] = GEP_PROB; // Uniform probability for GAPs
			else if (sequenceCode == seqCode) p[seqCode] = 1.0;
			else p[seqCode] = 0.0;
			return p[seqCode];
		}

		//---
		// Non-leaf node
		//---

		// Likelihood from the left sub-tree
		RealMatrix P = tmatrix.matrix(distanceLeft);
		double pleft = 0;
		if (left != null) {
			for (int sc = 0; sc < p.length; sc++) {
				double lleft = ((LikelihoodTree) left).likelihood(tmatrix, sc);
				double pij = P.getEntry(seqCode, sc);
				pleft += lleft * pij;
			}
		} else pleft = 1.0; // No node

		// Likelihood from the right sub-tree
		P = tmatrix.matrix(distanceRight);
		double pright = 0;
		if (right != null) {
			for (int sc = 0; sc < p.length; sc++) {
				double lright = ((LikelihoodTree) right).likelihood(tmatrix, sc);
				double pij = P.getEntry(seqCode, sc);
				pright += lright * pij;
			}
		} else pright = 1.0; // No node

		// Set code
		p[seqCode] = pleft * pright;

		return p[seqCode];
	}

	@Override
	protected PhylogeneticTree newNode(PhylogeneticTree parent, String phyloStr) {
		return new LikelihoodTree(parent, phyloStr);
	}

	@Override
	protected void resetNode() {
		if (p == null) p = new double[GprSeq.AMINO_ACIDS.length];
		Arrays.fill(p, Double.NaN);
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
