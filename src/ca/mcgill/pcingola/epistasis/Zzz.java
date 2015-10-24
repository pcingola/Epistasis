package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.coEvolutionMetrics.AaSimilarityMatrix;
import ca.mcgill.pcingola.epistasis.coEvolutionMetrics.McBasc;

public class Zzz {

	public static void main(String[] args) {
		String dir = Gpr.HOME + "/snpEff/epistasis_new";
		String matrix = dir + "/McLachlan_matrix.txt";

		AaSimilarityMatrix simMatrix = new AaSimilarityMatrix(matrix);
		System.out.println(simMatrix);

		String coli = "AGH";
		String colj = "AGH";
		McBasc mcBasc = new McBasc(simMatrix, coli, colj);
		System.out.println("Score: " + mcBasc.score());
	}

}
