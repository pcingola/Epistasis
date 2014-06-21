package ca.mcgill.pcingola.epistasis.entropy;

import gnu.trove.map.hash.TLongShortHashMap;
import gnu.trove.procedure.TLongShortProcedure;

import java.util.Arrays;

import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.MsaSimilarity;

public class Entropy {

	public static boolean debug = false;
	public static final double LOG_2 = Math.log(2.0);

	/**
	 * Conditional entropy
	 */
	public static double condEntropy(byte codei[], byte codej[]) {
		//---
		// Initialize counters
		//---
		int count = 0;
		int numAligns = codei.length;
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
			byte basei = codei[i];
			byte basej = codej[i];
			if ((basei == MsaSimilarity.ALIGN_GAP) || (basej == MsaSimilarity.ALIGN_GAP)) continue;

			// Count
			countI[basei]++;
			countJ[basej]++;
			countIJ[basei][basej]++;
			if (debug) System.out.println(GprSeq.code2aa(basei) + " / " + GprSeq.code2aa(basej) + "\t" + countIJ[basei][basej]);
			count++;
		}

		double hcondXY = 0.0;
		for (int i = 0; i < aaLen; i++) {
			if (countI[i] <= 0) continue;

			for (int j = 0; j < aaLen; j++) {
				if (countJ[j] <= 0) continue;

				if (countIJ[i][j] > 0) {
					double pij = (countIJ[i][j]) / ((double) count);
					double pi = (countI[i]) / ((double) count);
					double pj = (countJ[j]) / ((double) count);

					hcondXY += pij * Math.log(pj / pij) / LOG_2;

					if (debug) System.out.println(GprSeq.code2aa(codei[i]) + " / " + GprSeq.code2aa(codej[j]) + "\tpi: " + pi + "\tpj: " + pj + "\tpij: " + pij + "\th(x|y): " + hcondXY);
				}
			}
		}

		// Results
		return hcondXY;
	}

	/**
	 * Conditional entropy between AA sequences
	 */
	public static double condEntropy(String sequencei, String sequencej) {
		if (sequencei.length() != sequencej.length()) throw new RuntimeException("Lengths do not match!");

		// Convert string to byte codes
		int numAligns = sequencei.length();
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(sequencei.charAt(i));
			codej[i] = GprSeq.aa2Code(sequencej.charAt(i));
		}

		return condEntropy(codei, codej);
	}

	/**
	 * Entropy of a single sequence: H(X)
	 */
	public static double entropy(byte codei[]) {
		//---
		// Initialize counters
		//---
		int count = 0;
		int numAligns = codei.length;
		int aaLen = GprSeq.AA_TO_CODE.length;
		short countI[] = new short[aaLen];
		Arrays.fill(countI, (short) 0);

		//---
		// Count matching bases
		//---
		for (int i = 0; i < numAligns; i++) {
			byte basei = codei[i];
			if (basei == MsaSimilarity.ALIGN_GAP) continue;

			// Count
			countI[basei]++;
			count++;
		}

		double h = 0;
		for (int i = 0; i < aaLen; i++) {
			if (countI[i] <= 0) continue;
			if (countI[i] > 0) {
				double pi = (countI[i]) / ((double) count);
				h -= pi * Math.log(pi) / LOG_2;
			}
		}

		return h;
	}

	/**
	 * Entropy of two sequences: H(X,Y)
	 */
	public static double entropy(byte codei[], byte codej[]) {
		//---
		// Initialize counters
		//---
		int count = 0;
		int numAligns = codei.length;
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
			byte basei = codei[i];
			byte basej = codej[i];
			if ((basei == MsaSimilarity.ALIGN_GAP) || (basej == MsaSimilarity.ALIGN_GAP)) continue;

			// Count
			countI[basei]++;
			countJ[basej]++;
			countIJ[basei][basej]++;
			if (debug) System.out.println(GprSeq.code2aa(basei) + " / " + GprSeq.code2aa(basej) + "\t" + countIJ[basei][basej]);
			count++;
		}

		double hxy = 0;
		for (int i = 0; i < aaLen; i++) {
			if (countI[i] <= 0) continue;

			for (int j = 0; j < aaLen; j++) {
				if (countJ[j] <= 0) continue;

				if (countIJ[i][j] > 0) {
					double pij = (countIJ[i][j]) / ((double) count);
					double pi = (countI[i]) / ((double) count);
					double pj = (countJ[j]) / ((double) count);

					hxy -= pij * Math.log(pij) / LOG_2;

					if (debug) System.out.println(GprSeq.code2aa(codei[i]) + " / " + GprSeq.code2aa(codej[j]) + "\tpi: " + pi + "\tpj: " + pj + "\tpij: " + pij + "\thxy: " + hxy);
				}
			}
		}

		return hxy;
	}

	/**
	 * Entropy of a single sequence: H(X)
	 */
	public static double entropy(String coli) {
		// Convert string to byte codes
		int numAligns = coli.length();
		byte codei[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++)
			codei[i] = GprSeq.aa2Code(coli.charAt(i));

		return entropy(codei);
	}

	/**
	 * Entropy of two sequences: H(X,Y)
	 */
	public static double entropy(String coli, String colj) {
		if (coli.length() != colj.length()) throw new RuntimeException("Lengths do not match!");

		// Convert string to byte codes
		int numAligns = coli.length();
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(coli.charAt(i));
			codej[i] = GprSeq.aa2Code(colj.charAt(i));
		}

		return entropy(codei, codej);
	}

	/**
	 * Mutual Information
	 */
	public static double mutualInformation(byte codei[], byte codej[]) {
		//---
		// Initialize counters
		//---
		int count = 0;
		int numAligns = codei.length;
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
			byte basei = codei[i];
			byte basej = codej[i];
			if ((basei == MsaSimilarity.ALIGN_GAP) || (basej == MsaSimilarity.ALIGN_GAP)) continue;

			// Count
			countI[basei]++;
			countJ[basej]++;
			countIJ[basei][basej]++;
			count++;
		}

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
		return mutInf;
	}

	/**
	 * Mutual Information
	 */
	public static double mutualInformation(String coli, String colj) {
		if (coli.length() != colj.length()) throw new RuntimeException("Lengths do not match!");

		// Convert string to byte codes
		int numAligns = coli.length();
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(coli.charAt(i));
			codej[i] = GprSeq.aa2Code(colj.charAt(i));
		}

		return mutualInformation(codei, codej);
	}

	/**
	 * Variation of information = H(X|Y) + H(Y|X) = H(X,y) - I(X,Y)
	 * It can be used as a distance metric
	 */
	public static double variationOfInformation(byte codei[], byte codej[]) {
		//---
		// Initialize counters
		//---
		int count = 0;
		int numAligns = codei.length;
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
			byte basei = codei[i];
			byte basej = codej[i];
			if ((basei == MsaSimilarity.ALIGN_GAP) || (basej == MsaSimilarity.ALIGN_GAP)) continue;

			// Count
			countI[basei]++;
			countJ[basej]++;
			countIJ[basei][basej]++;
			count++;
		}

		double mutInf = 0.0;
		double hxy = 0;
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
		return hxy - mutInf;
	}

	public static double variationOfInformation(String coli, String colj) {
		if (coli.length() != colj.length()) throw new RuntimeException("Lengths do not match!");

		// Convert string to byte codes
		int numAligns = coli.length();
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(coli.charAt(i));
			codej[i] = GprSeq.aa2Code(colj.charAt(i));
		}

		return variationOfInformation(codei, codej);
	}

	final short ONE = 1;

	int count = 0;
	double mi;
	TLongShortHashMap countIJ = new TLongShortHashMap();
	TLongShortHashMap countI = new TLongShortHashMap();
	TLongShortHashMap countJ = new TLongShortHashMap();

	public Entropy() {
	}

	public long getCount() {
		return count;
	}

	public void inc(long basesI, long basesJ) {
		long basesIJ = (basesI << rot) | basesJ;

		countI.adjustOrPutValue(basesI, ONE, ONE);
		countJ.adjustOrPutValue(basesJ, ONE, ONE);
		countIJ.adjustOrPutValue(basesIJ, ONE, ONE);
		count++;
	}

	public double mi() {
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
