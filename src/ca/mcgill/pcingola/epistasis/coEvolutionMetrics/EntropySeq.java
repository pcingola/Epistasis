package ca.mcgill.pcingola.epistasis.coEvolutionMetrics;

import java.util.Arrays;

import ca.mcgill.mcb.pcingola.stats.BooleanMutable;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.msa.MsaSimilarity;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TObjectIntProcedure;

/**
 * Entropy and other functions for sequences
 *
 * @author pcingola
 */
public class EntropySeq {

	public enum InformationFunction {
		HXY, HCONDXY, MI, VARINF
	}

	public static final String SEPARATOR = "-";
	public static boolean debug = false;
	public static final double LOG_2 = Math.log(2.0);;

	final short ONE = 1;

	int count = 0;
	double mi, varInf, hxy, hx, hy, hcondXY, hcondYX;
	TObjectIntHashMap<String> countI = new TObjectIntHashMap<>();
	TObjectIntHashMap<String> countJ = new TObjectIntHashMap<>();
	TObjectIntHashMap<String> countIJ = new TObjectIntHashMap<>();

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
	 * Ratio: number of AA equal the the first / number of non-gap
	 */
	public static double conservation(String seq) {
		if (seq == null || seq.isEmpty()) return 0;

		char fisrt = seq.charAt(0);
		int countEq = 0, count = 0;
		for (int i = 1; i < seq.length(); i++) {
			if (seq.charAt(i) != '-') count++;
			if (seq.charAt(i) == fisrt) countEq++;

		}
		return ((double) countEq) / ((double) count);
	}

	/**
	 * Naive correlation
	 */
	public static double correlation(byte codei[], byte codej[]) {
		int count = 0, sum = 0;
		int len = codei.length;

		if (Math.random() < 2) throw new RuntimeException("THIS CALCULATION IS WRONG!");

		// Count matching bases
		for (int i = 0; i < len; i++) {
			byte basei = codei[i];
			byte basej = codej[i];

			// Do not take gaps into account
			if ((basei < 0) || (basej < 0)) continue;

			count++;
			if (basei == basej) sum++;
		}

		return ((double) sum) / ((double) count);
	}

	/**
	 * Naive correlation between two sequences
	 */
	public static double correlation(String seqi, String seqj) {
		// Convert string to byte codes
		int numAligns = seqi.length();
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(seqi.charAt(i));
			codej[i] = GprSeq.aa2Code(seqj.charAt(i));
		}

		return correlation(codei, codej);
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
		double varInf = hxy - mutInf;
		return varInf >= 0 ? varInf : 0.0;
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

	public EntropySeq() {
	}

	/**
	 * Calculate
	 */
	protected void calc() {
		hx = hy = hxy = hcondXY = hcondYX = mi = varInf = 0.0;
		BooleanMutable firstIteration = new BooleanMutable(true);
		if (count <= 0) return;

		countI.forEachEntry(new TObjectIntProcedure<String>() {

			@Override
			public boolean execute(String basesI, int countBasesI) {
				if (countBasesI <= 0) return true; // Nothing to do
				double pi = ((double) countBasesI) / (count);

				countJ.forEachEntry(new TObjectIntProcedure<String>() {

					@Override
					public boolean execute(String basesJ, int countBasesJ) {
						if (countBasesJ == 0) return true; // Nothing to do
						double pj = ((double) countBasesJ) / (count);
						if (firstIteration.is()) hy -= pj * Math.log(pj) / LOG_2;

						String basesIJ = basesI + SEPARATOR + basesJ;
						int countBasesIJ = countIJ.get(basesIJ);
						if (countBasesIJ == 0) return true; // Nothing to do

						double pij = ((double) countBasesIJ) / (count);
						mi += pij * Math.log(pij / (pi * pj)) / LOG_2;
						hxy -= pij * Math.log(pij) / LOG_2;
						hcondXY += pij * Math.log(pj / pij) / LOG_2;
						hcondYX += pij * Math.log(pi / pij) / LOG_2;

						return true;
					}
				});

				hx -= pi * Math.log(pi) / LOG_2;
				firstIteration.setFalse();

				return true;
			}
		});

		varInf = hxy - mi;
		if (varInf < 0) varInf = 0;
	}

	/**
	 * Calculate using sequences 'seqsi' and 'seqsj'
	 */
	public void calc(String[] seqsi, String[] seqsj) {
		// String length (all are assumed to be equal
		int len = 0;
		for (int i = 0; i < seqsi.length; i++)
			if (seqsi[i] != null) {
				len = seqsi[i].length();
				break;
			}

		// Count each alignment
		for (int i = 0; i < len; i++) {
			String basesi = seq(seqsi, i);
			String basesj = seq(seqsj, i);
			inc(basesi, basesj);
		}

		calc();
	}

	public long getCount() {
		return count;
	}

	public double getHcondXY() {
		return hcondXY;
	}

	public double getHcondYX() {
		return hcondYX;
	}

	public double getHx() {
		return hx;
	}

	public double getHxy() {
		return hxy;
	}

	public double getHy() {
		return hy;
	}

	public double getMi() {
		return mi;
	}

	public double getVarInf() {
		return varInf;
	}

	public void inc(String basesI, String basesJ) {
		countI.adjustOrPutValue(basesI, ONE, ONE);
		countJ.adjustOrPutValue(basesJ, ONE, ONE);
		countIJ.adjustOrPutValue(basesI + SEPARATOR + basesJ, ONE, ONE);
		count++;
	}

	/**
	 * Create a string from bases at position 'idx'
	 */
	String seq(String seqs[], int idx) {
		char bases[] = new char[seqs.length];
		for (int i = 0; i < seqs.length; i++)
			bases[i] = seqs[i] != null ? seqs[i].charAt(idx) : ' ';

		return new String(bases);
	}

}
