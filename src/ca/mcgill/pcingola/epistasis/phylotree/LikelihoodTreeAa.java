package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.Set;

import org.apache.commons.math3.linear.RealMatrix;

import ca.mcgill.mcb.pcingola.stats.Counter;
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

	public static final double GAP_PROB = 1.0;
	double p[];

	/**
	 * Create root node
	 */
	public LikelihoodTreeAa() {
		super();
	}

	public LikelihoodTreeAa(PhylogeneticTree parent, String phyloStr, Counter ids) {
		super(parent, phyloStr, ids);
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
		uniformCode();

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
		if (uniformCode == -1) return GAP_PROB;
		if (!Double.isNaN(p[aaCode])) return p[aaCode];

		//---
		// Leaf node?
		//---
		if (isLeaf()) {
			if (isGap()) p[aaCode] = GAP_PROB;// GAP probability

			// Probability is 1 for that sequence, 0 for others
			Arrays.fill(p, 0.0);
			p[sequenceCode] = 1.0;

			return p[aaCode];
		}

		//---
		// Non-leaf node
		//---
		if (uniformCode >= 0) {
			Double punif = lcache.get(this, aaCode);
			if (punif != null) return punif;
		}

		// Likelihood from the left sub-tree
		RealMatrix P = tmatrix.matrix(distanceLeft);
		double pleft = 0.0;
		if (left != null) {
			if (left.isLeaf()) {
				// Leaf node likelihoods are zero for all elements except 'sequenceCode' (which has a likelihood of 1.0)
				if (left.isGap()) pleft = 1.0;
				else pleft = P.getEntry(aaCode, left.sequenceCode);
			} else {
				// Sum likelihoods over all possible 'aa'
				for (int aa2 = 0; aa2 < p.length; aa2++) {
					double lleft = ((LikelihoodTreeAa) left).likelihood(tmatrix, aa2);
					double pij = P.getEntry(aaCode, aa2);
					pleft += lleft * pij;
				}
			}
		} else pleft = 1.0;

		// Likelihood from the right sub-tree
		P = tmatrix.matrix(distanceRight);
		double pright = 0.0;
		if (right != null) {
			if (right.isLeaf()) {
				// Leaf node likelihoods are zero for all elements except 'sequenceCode' (which has a likelihood of 1.0)
				if (right.isGap()) pright = 1.0;
				else pright = P.getEntry(aaCode, right.sequenceCode);
			} else {
				// Sum likelihoods over all possible 'aa'
				for (int aa2 = 0; aa2 < p.length; aa2++) {
					double lright = ((LikelihoodTreeAa) right).likelihood(tmatrix, aa2);
					double pij = P.getEntry(aaCode, aa2);
					pright += lright * pij;
				}
			}
		} else pright = 1.0;

		// Set code
		p[aaCode] = pleft * pright;

		// Update cache
		if (uniformCode >= 0) lcache.set(this, aaCode, p[aaCode]);

		return p[aaCode];
	}

	@Override
	protected PhylogeneticTree newNode(PhylogeneticTree parent, String phyloStr, Counter ids) {
		return new LikelihoodTreeAa(parent, phyloStr, ids);
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
		if (p.length != size) throw new RuntimeException("Calculating for different size not allowed!");

		Arrays.fill(p, Double.NaN);

		if (left != null) left.resetNode(size);
		if (right != null) right.resetNode(size);
	}

	/**
	 * Pre-calculate matrix exponentials
	 */
	public void times(Set<Double> times) {
		if (isLeaf()) return;

		times.add(distanceLeft);
		((LikelihoodTreeAa) left).times(times);

		times.add(distanceRight);
		((LikelihoodTreeAa) right).times(times);
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
