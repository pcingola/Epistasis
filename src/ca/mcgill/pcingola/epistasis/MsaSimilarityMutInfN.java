package ca.mcgill.pcingola.epistasis;

import gnu.trove.map.hash.TLongShortHashMap;
import gnu.trove.procedure.TLongShortProcedure;

/**
 * Calculate Mutual Information
 * 
 * @author pcingola
 */
class MICalc {

	public static final double LOG_2 = Math.log(2.0);
	final short ONE = 1;

	int count = 0;
	int num, rot;
	double mi;
	TLongShortHashMap countIJ = new TLongShortHashMap();
	TLongShortHashMap countI = new TLongShortHashMap();
	TLongShortHashMap countJ = new TLongShortHashMap();

	public MICalc(int num) {
		this.num = num;
		rot = MultipleSequenceAlignment.ROTATE_BITS * num;
	}

	void inc(long basesI, long basesJ) {
		long basesIJ = (basesI << rot) | basesJ;

		countI.adjustOrPutValue(basesI, ONE, ONE);
		countJ.adjustOrPutValue(basesJ, ONE, ONE);
		countIJ.adjustOrPutValue(basesIJ, ONE, ONE);
		count++;
	}

	double mi() {
		if (count <= 0) return 0.0;
		mi = 0.0;

		countI.forEachEntry(new TLongShortProcedure() {

			@Override
			public boolean execute(long basesI, short countBasesI) {
				if (countBasesI <= 0) return true; // Nothing to do

				countJ.forEachEntry(new TLongShortProcedure() {

					@Override
					public boolean execute(long basesJ, short countBasesJ) {
						if (countBasesJ == 0) return true; // Nothing to do

						long basesIJ = (basesI << rot) | basesJ;
						int countBasesIJ = countIJ.get(basesIJ);

						if (countBasesIJ == 0) return true; // Nothing to do
						// System.err.println("basesI: " + basesI + "\tcountBasesI: " + countBasesI + "\tbasesJ: " + basesJ + "\tcountBasesJ: " + countBasesJ + "\tcountBasesIJ: " + countBasesIJ);

						double pij = ((double) countBasesIJ) / (count);
						double pi = ((double) countBasesI) / (count);
						double pj = ((double) countBasesJ) / (count);
						mi += pij * Math.log(pij / (pi * pj)) / LOG_2;

						return true;
					}
				});

				return true;
			}
		});

		return mi;
	}
}

/**
 * Implement a 'similarity' by mutual information using N bases around the target
 * 
 * @author pcingola
 */
public class MsaSimilarityMutInfN extends MsaSimilarity {

	public static final double LOG_2 = Math.log(2.0);

	public MsaSimilarityMutInfN(MultipleSequenceAlignmentSet msas, int numBases) {
		super(msas);
		maxScore = 5.0;
		this.numBases = numBases;
	}

	/**
	 * Measure similarity: Correlation between two loci
	 */
	@Override
	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		// Sanity check
		if ((posi + numBases) >= msai.getSeqLen()) return Double.NaN;
		if ((posj + numBases) >= msaj.getSeqLen()) return Double.NaN;

		int numAligns = msas.getNumAligns();

		//---
		// Initialize counters
		//---
		MICalc miCalc = new MICalc(numBases);

		//---
		// Count matching bases
		//---
		for (int i = 0; i < numAligns; i++) {
			long basesI = msai.getCodeLong(i, posi, numBases);
			if (basesI == ALIGN_GAP) continue;

			long basesJ = msaj.getCodeInt(i, posj, numBases);
			if (basesJ == ALIGN_GAP) continue;

			// Count
			miCalc.inc(basesI, basesJ);
		}

		// Only a few organisms align? Filter out
		if (miCalc.count < MIN_COUNT_THRESHOLD) return Double.NaN;

		double mutInf = miCalc.mi();

		// Results
		incScore(mutInf);
		if (mutInf > threshold) return mutInf;
		return Double.NaN;
	}
}
