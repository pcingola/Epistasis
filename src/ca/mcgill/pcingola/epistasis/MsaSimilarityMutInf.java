package ca.mcgill.pcingola.epistasis;

import java.util.Arrays;

import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Implement a 'similarity' by mutual information
 *
 * @author pcingola
 */
public class MsaSimilarityMutInf extends MsaSimilarity {

	public static final double LOG_2 = Math.log(2.0);

	public MsaSimilarityMutInf(MultipleSequenceAlignmentSet msas) {
		super(msas);
		double n = GprSeq.AMINO_ACIDS.length;
		double p = 1.0 / n;
		maxScore = -Math.log(p) / Math.log(2.0); // Maximum possible entropy
	}

	/**
	 * Measure similarity: Mutual information
	 */
	@Override
	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		int count = 0;
		int numAligns = msas.getNumAligns();

		//---
		// Initialize counters
		//---
		int aaLen = GprSeq.AA_TO_CODE.length;
		short countIJ[][] = new short[aaLen][aaLen];
		short countI[] = new short[aaLen];
		short countJ[] = new short[aaLen];

		for (int i = 0; i < aaLen; i++)
			Arrays.fill(countIJ[i], (short) 0);
		Arrays.fill(countI, (short) 0);
		Arrays.fill(countJ, (short) 0);

		//---
		// Count matching bases
		//---
		for (int i = 0; i < numAligns; i++) {
			byte basei = msai.getCode(i, posi);
			byte basej = msaj.getCode(i, posj);
			if ((basei < 0) || (basej < 0)) continue;

			// Count
			countI[basei]++;
			countJ[basej]++;
			countIJ[basei][basej]++;
			count++;
		}

		// Only a few organisms align? Filter out
		if (count < minCount) return Double.NaN;

		double mutInf = 0.0;
		for (int i = 0; i < aaLen; i++) {
			if (countI[i] <= 0) continue;

			for (int j = 0; j < aaLen; j++) {
				if (countJ[j] <= 0) continue;

				if (countIJ[i][j] > 0) {
					double pij = ((double) countIJ[i][j]) / ((double) count);
					double pi = ((double) countI[i]) / ((double) count);
					double pj = ((double) countJ[j]) / ((double) count);

					mutInf += pij * Math.log(pij / (pi * pj)) / LOG_2;
				}
			}
		}

		// Results
		incScore(mutInf);
		if (mutInf >= threshold) return mutInf;
		return Double.NaN;
	}

}
