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
public class TestCaseZzz extends TestCase {

	public static boolean debug = true;
	public static final double EPSILON = 1E-6;

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
