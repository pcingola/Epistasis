package ca.mcgill.pcingola.epistasis.phylotree;

import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.DistanceResults;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignmentSet;

/**
 * Estimate  transition matrix of AA Pairs
 */
public class EstimateTransitionMatrixPairs extends EstimateTransitionMatrix {

	DistanceResults aaContacts;

	public EstimateTransitionMatrixPairs(LikelihoodTree tree, MultipleSequenceAlignmentSet msas, DistanceResults aaContacts) {
		super(tree, msas);
		this.aaContacts = aaContacts;
		N = GprSeq.AMINO_ACIDS.length * GprSeq.AMINO_ACIDS.length;
		verbose = true;
	}

	/**
	 * Count occurrences
	 */
	@Override
	protected int[] countAa() {
		return msas.countAaPairs(aaContacts);
	}

	/**
	 *  Count all transitions
	 */
	@Override
	protected int[][] countTransitions(int seqNum1, int seqNum2) {
		return msas.countTransitionsPairs(seqNum1, seqNum2, aaContacts);
	}

	@Override
	protected void initNames() {
		// Column and row names
		names = new String[GprSeq.AMINO_ACIDS.length * GprSeq.AMINO_ACIDS.length];
		for (int i = 0, h = 0; i < GprSeq.AMINO_ACIDS.length; i++)
			for (int j = 0; j < GprSeq.AMINO_ACIDS.length; j++, h++)
				names[h] = GprSeq.code2aa((byte) i) + "_" + GprSeq.code2aa((byte) j);
	}

}
