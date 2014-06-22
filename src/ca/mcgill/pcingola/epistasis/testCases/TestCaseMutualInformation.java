package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.pcingola.epistasis.MsaSimilarityMutInf;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.entropy.EntropySeq;

/**
 * Implement a 'similarity' by mutual information and conditional entropy
 *
 * @author pcingola
 */
public class TestCaseMutualInformation extends TestCase {

	public static final double EPSILON = 1E-6;
	public static final double LOG_2 = Math.log(2.0);
	public static boolean debug = false;

	/**
	 * Measure similarity: Correlation between two loci
	 */
	public static double z(String coli, String colj) {
		double hxy = EntropySeq.entropy(coli, colj);
		double hcondXY = EntropySeq.condEntropy(coli, colj);
		double hcondYX = EntropySeq.condEntropy(colj, coli);

		System.out.println("h(x,y): " + hxy + "\th(x|y): " + hcondXY + "\th(y|x): " + hcondYX);
		// Results
		return hxy + (hcondXY + hcondYX);
	}

	public void checkEntropy(String seqi, String seqj) {
		// Create "multiple columns" (but just use only one)
		String colsi[] = { seqi };
		String colsj[] = { seqj };
		checkEntropy(colsi, colsj);
	}

	public void checkEntropy(String colsi[], String colsj[]) {
		EntropySeq entropy = new EntropySeq();
		entropy.calc(colsi, colsj);

		System.out.println(colsi[0] + "\n" + colsj[0] //
				+ "\n\tH(X)        = " + EntropySeq.entropy(colsi[0]) + "\t" + entropy.getHx() //
				+ "\n\tH(Y)        = " + EntropySeq.entropy(colsj[0]) + "\t" + entropy.getHy() //
				+ "\n\tH(X,Y)      = " + EntropySeq.entropy(colsi[0], colsj[0]) + "\t" + entropy.getHxy() //
				+ "\n\tH(X|Y)      = " + EntropySeq.condEntropy(colsi[0], colsj[0]) + "\t" + entropy.getHcondXY() //
				+ "\n\tH(Y|X)      = " + EntropySeq.condEntropy(colsj[0], colsi[0]) + "\t" + entropy.getHcondYX() //
				+ "\n\tMI          = " + EntropySeq.mutualInformation(colsi[0], colsj[0]) + "\t" + entropy.getMi() //
				+ "\n\tVarInf(X,Y) = " + EntropySeq.variationOfInformation(colsj[0], colsi[0]) + "\t" + entropy.getVarInf() //
				);

		Assert.assertEquals(EntropySeq.mutualInformation(colsi[0], colsj[0]), entropy.getMi(), 1E-6);
		Assert.assertEquals(EntropySeq.entropy(colsi[0]), entropy.getHx(), 1E-6);
		Assert.assertEquals(EntropySeq.entropy(colsj[0]), entropy.getHy(), 1E-6);
		Assert.assertEquals(EntropySeq.entropy(colsi[0], colsj[0]), entropy.getHxy(), 1E-6);
		Assert.assertEquals(EntropySeq.condEntropy(colsi[0], colsj[0]), entropy.getHcondXY(), 1E-6);
		Assert.assertEquals(EntropySeq.condEntropy(colsj[0], colsi[0]), entropy.getHcondYX(), 1E-6);
		Assert.assertEquals(EntropySeq.variationOfInformation(colsj[0], colsi[0]), entropy.getVarInf(), 1E-6);
	}

	public void checkEntropyMultiple(String seqi, String seqj) {
		// Create "multiple columns" (but just use only one)
		String colsi[] = { seqi, seqi, seqi };
		String colsj[] = { seqj, seqj, seqj };
		checkEntropy(colsi, colsj);
	}

	public void test_01() {
		String coli = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";
		String colj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL";

		double mi = EntropySeq.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(1.0, mi, 1E-6);
	}

	public void test_02() {
		String coli = "AAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDCCCCCCCCCCCCCCCCCCCC";
		String colj = "KKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLPPPPPPPPPPPPPPPPPPPP";

		double mi = EntropySeq.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(1.584962500721156, mi, 1E-6);
	}

	public void test_03() {
		String coli = "ARNDCEQGHILKMFPSTWYV";
		String colj = "RNDCEQGHILKMFPSTWYVA";

		double mi = EntropySeq.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(4.321928094887363, mi, 1E-6);
	}

	public void test_04() {
		String coli = "ARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYV";
		String colj = "RNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVA";

		double mi = EntropySeq.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(4.321928094887363, mi, 1E-6);
	}

	public void test_05() {
		String coli = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDAA";
		String colj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLKK";

		double mi = EntropySeq.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(1.0, mi, 1E-6);
	}

	public void test_06() {
		MultipleSequenceAlignmentSet msas = new MultipleSequenceAlignmentSet("test/msa_1.fa", 100);
		msas.load();

		MsaSimilarityMutInf simMi = new MsaSimilarityMutInf(msas);
		MultipleSequenceAlignment msa = msas.getMsas().get(0);

		// Calculate MI for every pair
		for (int i = 0; i < msa.length(); i++) {
			for (int j = i + 1; j < msa.length(); j++) {
				String seqi = msa.getColumnString(i);
				String seqj = msa.getColumnString(j);

				// Calculate MI both ways and compare
				double mi = EntropySeq.mutualInformation(seqi, seqj);
				double mi2 = simMi.calc(msa, msa, i, j);

				if (debug) System.out.println("sequences: " + i + " , " + j + "\tmi: " + mi + "\tmi2: " + mi2 + "\n\t" + seqi + "\n\t" + seqj + "\n");
				Assert.assertEquals(mi, mi2, EPSILON);
			}
		}
	}

	public void test_11() {
		String seqi = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";
		String seqj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL";
		checkEntropy(seqi, seqj);
	}

	public void test_12() {
		String seqi = "AAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDCCCCCCCCCCCCCCCCCCCC";
		String seqj = "KKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLPPPPPPPPPPPPPPPPPPPP";
		checkEntropy(seqi, seqj);
	}

	public void test_13() {
		String seqi = "ARNDCEQGHILKMFPSTWYV";
		String seqj = "RNDCEQGHILKMFPSTWYVA";
		checkEntropy(seqi, seqj);
	}

	public void test_14() {
		String seqi = "ARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYV";
		String seqj = "RNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVA";
		checkEntropy(seqi, seqj);
	}

	public void test_15() {
		String seqi = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDAA";
		String seqj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLKK";
		checkEntropy(seqi, seqj);
	}

	public void test_20() {
		String seqi = "ARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYV";
		String seqj = "RNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVA";
		checkEntropyMultiple(seqi, seqj);
	}

}
