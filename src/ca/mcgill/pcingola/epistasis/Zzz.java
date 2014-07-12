package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrix;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrixMarkov;

/**
 * Test
 *
 * @author pcingola
 */
public class Zzz {

	public static final boolean debug = false;

	public static void main(String[] args) {
		TransitionMatrix qhat = TransitionMatrixMarkov.load(Gpr.HOME + "/snpEff/epistasis/Qhat.txt");
		System.out.println(qhat.toStringNice());
	}
}
