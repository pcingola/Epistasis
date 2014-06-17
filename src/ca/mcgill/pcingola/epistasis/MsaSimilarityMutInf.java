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

	/**
	 * Measure similarity: Mutual information (do not return NaN, use 1.0 instead)
	 */
	public static double miNoNan(String seqi, String seqj) {
		double mi = mi(seqi, seqj);
		return Double.isNaN(mi) ? 0.0 : mi;
	}

	/**
	 * Measure similarity: Mutual information
	 */
	public static double mi(String seqi, String seqj) {
		int count = 0;

		if (seqi.length() != seqj.length()) throw new RuntimeException("Sequences length differ\n\t" + seqi + "\n\t" + seqj);

		//---
		// Initialize counters
		//---
		int numAa = GprSeq.AA_TO_CODE.length;
		short countIJ[][] = new short[numAa][numAa];
		short countI[] = new short[numAa];
		short countJ[] = new short[numAa];
		byte basesI[] = new byte[seqi.length()];
		byte basesJ[] = new byte[seqj.length()];

		for (int i = 0; i < numAa; i++)
			Arrays.fill(countIJ[i], (short) 0);

		int aaLen = seqi.length();
		for (int i = 0; i < basesI.length; i++) {
			basesI[i] = GprSeq.aa2Code(seqi.charAt(i));
			basesJ[i] = GprSeq.aa2Code(seqj.charAt(i));
		}

		Arrays.fill(countI, (short) 0);
		Arrays.fill(countJ, (short) 0);

		//---
		// Count matching bases
		//---
		for (int i = 0; i < aaLen; i++) {
			byte basei = basesI[i];
			byte basej = basesJ[i];
			if ((basei < 0) || (basej < 0)) continue;

			// Count
			countI[basei]++;
			countJ[basej]++;
			countIJ[basei][basej]++;
			count++;
		}

		// Only a few organisms align? Filter out
		if (count < MIN_COUNT_THRESHOLD) return Double.NaN;

		double mutInf = 0.0;
		for (int i = 0; i < numAa; i++) {
			if (countI[i] <= 0) continue;

			for (int j = 0; j < numAa; j++) {
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
		return mutInf;
	}

	protected double threshold = 1.5;

	public MsaSimilarityMutInf(MultipleSequenceAlignmentSet msas) {
		super(msas);
		maxScore = 4.0;
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
		if (count < MIN_COUNT_THRESHOLD) return Double.NaN;

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
		if (mutInf > threshold) return mutInf;
		return Double.NaN;
	}

}
