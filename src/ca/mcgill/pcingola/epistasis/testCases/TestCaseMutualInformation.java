package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.pcingola.epistasis.MsaSimilarityMutInf;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.entropy.Entropy;

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
		double hxy = Entropy.entropy(coli, colj);
		double hcondXY = Entropy.condEntropy(coli, colj);
		double hcondYX = Entropy.condEntropy(colj, coli);

		System.out.println("h(x,y): " + hxy + "\th(x|y): " + hcondXY + "\th(y|x): " + hcondYX);
		// Results
		return hxy + (hcondXY + hcondYX);
	}

	public void test_01() {
		String coli = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";
		String colj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL";

		double mi = Entropy.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(1.0, mi, 1E-6);
	}

	public void test_02() {
		String coli = "AAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDCCCCCCCCCCCCCCCCCCCC";
		String colj = "KKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLPPPPPPPPPPPPPPPPPPPP";

		double mi = Entropy.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(1.584962500721156, mi, 1E-6);
	}

	public void test_03() {
		String coli = "ARNDCEQGHILKMFPSTWYV";
		String colj = "RNDCEQGHILKMFPSTWYVA";

		double mi = Entropy.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(4.321928094887363, mi, 1E-6);
	}

	public void test_04() {
		String coli = "ARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYV";
		String colj = "RNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVA";

		double mi = Entropy.mutualInformation(coli, colj);
		double hsum = z(coli, colj);
		System.out.println("MI: " + mi + "\th_sum: " + hsum + "\n\t" + coli + "\n\t" + colj + "\n");
		Assert.assertEquals(4.321928094887363, mi, 1E-6);
	}

	public void test_05() {
		String coli = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDAA";
		String colj = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLKK";

		double mi = Entropy.mutualInformation(coli, colj);
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
				double mi = Entropy.mutualInformation(seqi, seqj);
				double mi2 = simMi.calc(msa, msa, i, j);

				if (debug) System.out.println("sequences: " + i + " , " + j + "\tmi: " + mi + "\tmi2: " + mi2 + "\n\t" + seqi + "\n\t" + seqj + "\n");
				Assert.assertEquals(mi, mi2, EPSILON);
			}
		}
	}

}
