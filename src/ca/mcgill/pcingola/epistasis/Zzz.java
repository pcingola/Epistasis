package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.likelihood.ParameterDistributionModel;

public class Zzz {

	public static void main(String[] args) {
		String file = Gpr.HOME + "/snpEff/epistasis/gwas/gwas.30.head.txt";
		String gwasTheta = Gpr.HOME + "/snpEff/epistasis/gwas/gwasCoeff_altNull.txt";

		ParameterDistributionModel pdmAlt = new ParameterDistributionModel(gwasTheta, "ALT");
		Gpr.debug(pdmAlt);

		ParameterDistributionModel pdmNull = new ParameterDistributionModel(gwasTheta, "NULL");
		Gpr.debug(pdmNull);

		//		Genome genome = new Genome("test");
		//		LineFileIterator lfi = new LineFileIterator(file);
		//		for (String line : lfi) {
		//			System.out.println(line);
		//			GwasResult gr = new GwasResult(genome, line);
		//			System.out.println("\t" + gr + "\n");
		//		}
	}
}
