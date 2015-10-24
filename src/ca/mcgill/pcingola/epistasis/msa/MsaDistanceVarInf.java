package ca.mcgill.pcingola.epistasis.msa;

import java.util.Arrays;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.coEvolutionMetrics.EntropySeq.InformationFunction;

/**
 * Implement a 'distance' by variation of information
 *
 * @author pcingola
 */
public class MsaDistanceVarInf extends MsaSimilarity {

	public MsaDistanceVarInf(MultipleSequenceAlignmentSet msas) {
		super(msas, InformationFunction.VARINF);
		double n = GprSeq.AMINO_ACIDS.length;
		double p = 1.0 / n;
		maxScore = 2.0 * -Math.log(p) / Math.log(2.0); // Twice the maximum entropy?
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
		double hxy = 0.0;
		for (int i = 0; i < aaLen; i++) {
			if (countI[i] <= 0) continue;

			for (int j = 0; j < aaLen; j++) {
				if (countJ[j] <= 0) continue;

				if (countIJ[i][j] > 0) {
					double pij = ((double) countIJ[i][j]) / ((double) count);
					double pi = ((double) countI[i]) / ((double) count);
					double pj = ((double) countJ[j]) / ((double) count);

					mutInf += pij * Math.log(pij / (pi * pj)) / LOG_2;
					hxy -= pij * Math.log(pij) / LOG_2;
				}
			}
		}

		// Results
		double varInf = hxy - mutInf;

		if (debug) {
			Gpr.debug("Zero!\th(x,y):" + hxy + "\tmi: " + mutInf //
					+ "\n\t" + msai.getId() + "[" + posi + "]:\t" + msai.getColumnString(posi) + "\t" + msai.isSkip(posi) //
					+ "\n\t" + msaj.getId() + "[" + posj + "]:\t" + msaj.getColumnString(posj) + "\t" + msaj.isSkip(posj) //
					+ "\n");
		}
		incScore(varInf);
		return varInf;
	}

}
