package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.pcingola.epistasis.entropy.Entropy;

/**
 * Implement a 'similarity' by mutual information and conditional entropy
 *
 * @author pcingola
 */
public class TestCaseZzz extends TestCase {

	public static boolean debug = true;
	public static final double EPSILON = 1E-6;

	public void checkEntropyMultiple(String seqi, String seqj) {
		// Create "multiple columns" (but just use only one)
		String colsi[] = { seqi, seqi, seqi };
		String colsj[] = { seqj, seqj, seqj };
		checkEntropy(colsi, colsj);
	}

	public void checkEntropy(String colsi[], String colsj[]) {
		Entropy entropy = new Entropy();
		entropy.calc(colsi, colsj);

		System.out.println(colsi[0] + "\n" + colsj[0] //
				+ "\n\tH(X)        = " + Entropy.entropy(colsi[0]) + "\t" + entropy.getHx() //
				+ "\n\tH(Y)        = " + Entropy.entropy(colsj[0]) + "\t" + entropy.getHy() //
				+ "\n\tH(X,Y)      = " + Entropy.entropy(colsi[0], colsj[0]) + "\t" + entropy.getHxy() //
				+ "\n\tH(X|Y)      = " + Entropy.condEntropy(colsi[0], colsj[0]) + "\t" + entropy.getHcondXY() //
				+ "\n\tH(Y|X)      = " + Entropy.condEntropy(colsj[0], colsi[0]) + "\t" + entropy.getHcondYX() //
				+ "\n\tMI          = " + Entropy.mutualInformation(colsi[0], colsj[0]) + "\t" + entropy.getMi() //
				+ "\n\tVarInf(X,Y) = " + Entropy.variationOfInformation(colsj[0], colsi[0]) + "\t" + entropy.getVarInf() //
		);

		Assert.assertEquals(Entropy.mutualInformation(colsi[0], colsj[0]), entropy.getMi(), 1E-6);
		Assert.assertEquals(Entropy.entropy(colsi[0]), entropy.getHx(), 1E-6);
		Assert.assertEquals(Entropy.entropy(colsj[0]), entropy.getHy(), 1E-6);
		Assert.assertEquals(Entropy.entropy(colsi[0], colsj[0]), entropy.getHxy(), 1E-6);
		Assert.assertEquals(Entropy.condEntropy(colsi[0], colsj[0]), entropy.getHcondXY(), 1E-6);
		Assert.assertEquals(Entropy.condEntropy(colsj[0], colsi[0]), entropy.getHcondYX(), 1E-6);
		Assert.assertEquals(Entropy.variationOfInformation(colsj[0], colsi[0]), entropy.getVarInf(), 1E-6);
	}

	public void test_20() {
		String seqi = "ARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYV";
		String seqj = "RNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVARNDCEQGHILKMFPSTWYVA";
		checkEntropyMultiple(seqi, seqj);
	}

}
