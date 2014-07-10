package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.phylotree.TransitionMatrix;

/**
 * Test
 *
 * @author pcingola
 */
public class Zzz {

	public static final boolean debug = false;

	public static void main(String[] args) {
		TransitionMatrix qhat = new TransitionMatrix(Gpr.HOME + "/snpEff/z.txt");
		System.out.println(qhat);
	}
}
