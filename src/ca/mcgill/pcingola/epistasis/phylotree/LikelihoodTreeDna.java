package ca.mcgill.pcingola.epistasis.phylotree;

import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * A phylogenetic tree that can be used to calculate likelihoods
 *
 * Ref: Felsestein algorithm
 * 		"Computational Molecular Evolution" Z. Yang, Section 4.2.2, page 102
 * 		"Biological Sequence Analysis", Durbin et. al., page 199
 *
 * @author pcingola
 */
public class LikelihoodTreeDna extends LikelihoodTreeAa {

	/**
	 * Create root node
	 */
	public LikelihoodTreeDna() {
		super();
	}

	public LikelihoodTreeDna(String name) {
		super(name);
	}

	public LikelihoodTreeDna(String name, PhylogeneticTree left, double distanceLeft, PhylogeneticTree right, double distanceRight) {
		super(name, left, distanceLeft, right, distanceRight);
	}

	@Override
	public char getSequence() {
		if (sequenceCode < 0) return ' ';
		return GprSeq.code2dna((byte) sequenceCode);
	}

	@Override
	public void setSequence(char dna) {
		sequenceCode = GprSeq.dna2Code(dna);
	}

}
