package ca.mcgill.pcingola.epistasis.testCases;

import java.util.Arrays;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.MsaSimilarity;

/**
 * Implement a 'similarity' by mutual information and conditional entropy
 *
 * @author pcingola
 */
public class TestCaseMutualInformation extends TestCase {

	public static final double LOG_2 = Math.log(2.0);
	public static boolean debug = false;

	/**
	 * Measure similarity: Correlation between two loci
	 */
	public static double condEntropy(String coli, String colj) {
		if (coli.length() != colj.length()) throw new RuntimeException("Lengths do not match!");

		int count = 0;
		int numAligns = coli.length();

		// Convert string to byte codes
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(coli.charAt(i));
			codej[i] = GprSeq.aa2Code(colj.charAt(i));
		}

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
					double pij = ((double) countIJ[i][j]) / ((double) count);
					double pi = ((double) countI[i]) / ((double) count);
					double pj = ((double) countJ[j]) / ((double) count);

					hcondXY += pij * Math.log(pj / pij) / LOG_2;

					if (debug) System.out.println(GprSeq.code2aa(codei[i]) + " / " + GprSeq.code2aa(codej[j]) + "\tpi: " + pi + "\tpj: " + pj + "\tpij: " + pij + "\th(x|y): " + hcondXY);
				}
			}
		}

		// Results
		return hcondXY;
	}

	/**
	 * Entropy
	 */
	public static double entropy(String coli, String colj) {

		if (coli.length() != colj.length()) throw new RuntimeException("Lengths do not match!");

		int count = 0;
		int numAligns = coli.length();

		// Convert string to byte codes
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(coli.charAt(i));
			codej[i] = GprSeq.aa2Code(colj.charAt(i));
		}

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
					double pij = ((double) countIJ[i][j]) / ((double) count);
					double pi = ((double) countI[i]) / ((double) count);
					double pj = ((double) countJ[j]) / ((double) count);

					hxy -= pij * Math.log(pij) / LOG_2;

					if (debug) System.out.println(GprSeq.code2aa(codei[i]) + " / " + GprSeq.code2aa(codej[j]) + "\tpi: " + pi + "\tpj: " + pj + "\tpij: " + pij + "\thxy: " + hxy);
				}
			}
		}

		return hxy;
	}

	/**
	 * Measure similarity: Mutual Information
	 */
	public static double mutualInformation(String coli, String colj) {

		if (coli.length() != colj.length()) throw new RuntimeException("Lengths do not match!");

		int count = 0;
		int numAligns = coli.length();

		// Convert string to byte codes
		byte codei[] = new byte[numAligns];
		byte codej[] = new byte[numAligns];
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(coli.charAt(i));
			codej[i] = GprSeq.aa2Code(colj.charAt(i));
		}

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
	 * Measure similarity: Correlation between two loci
	 */
	public static double z(String coli, String colj) {
		double hxy = entropy(coli, colj);
		double hcondXY = condEntropy(coli, colj);
		double hcondYX = condEntropy(colj, coli);

		System.out.println("h(x,y): " + hxy + "\th(x|y): " + hcondXY + "\th(y|x): " + hcondYX);
		// Results
		return hxy + (hcondXY + hcondYX);
	}

	public void test_01() {
		String coli = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";
		String colj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL";

		double mi = mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj);
		Assert.assertEquals(1.0, mi, 1E-6);
	}

	public void test_02() {
		String coli = "AAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDCCCCCCCCCCCCCCCCCCCC";
		String colj = "KKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLPPPPPPPPPPPPPPPPPPPP";

		double mi = mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj);
		Assert.assertEquals(1.584962500721156, mi, 1E-6);
	}

	public void test_03() {
		String coli = "ARNDCEQGHILKMFPSTWYV";
		String colj = "RNDCEQGHILKMFPSTWYVA";

		double mi = mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj);
		Assert.assertEquals(4.321928094887363, mi, 1E-6);
	}

	public void test_04() {
		String coli = "ARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYV";
		String colj = "RNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVA";

		double mi = mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj);
		Assert.assertEquals(4.321928094887363, mi, 1E-6);
	}

	public void test_05() {
		String coli = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDAA";
		String colj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLKK";

		double mi = mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj);
		Assert.assertEquals(1.0, mi, 1E-6);
	}

}
