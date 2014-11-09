package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.GenotypePos;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.pdb.PdbGenomeMsas;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTreeAa;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseZzz extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = true || debug;

	public void test_zzz() {
		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String phyloFileName = "data/hg19.100way.nh";
		String msasFile = Gpr.HOME + "/snpEff/epistasis/msas.best.fa"; // "data/msa_test.fa.gz";
		String pdbDir = ""; // Not used

		LikelihoodTreeAa tree = new LikelihoodTreeAa();
		tree.load(phyloFileName);

		MultipleSequenceAlignmentSet msas = new MultipleSequenceAlignmentSet(msasFile, tree);
		msas.load();
		msas.buildForest();

		PdbGenomeMsas pdbGenomeMsas = new PdbGenomeMsas(configFile, genome, pdbDir, msas);
		pdbGenomeMsas.setDebug(debug);
		pdbGenomeMsas.setTree(tree);
		pdbGenomeMsas.initialize();

		/**
		 * VCF entry:
		 * 6	31864410	rs149384831	C	T
		 * 
		 * Error:
		 *         ID                : 6:31864409_C/T
		 *         Marker            : Marker_6:31864383-31864601
		 *         msa.Id            : NM_006709_6:31864382-31864600
		 *         msa.aaIdx         : -1
		 */
		Chromosome chr = pdbGenomeMsas.getConfig().getGenome().getOrCreateChromosome("6");
		GenotypePos gp = new GenotypePos(chr, 31864409, "6:31864409_C/T");
		gp.mapGenomic2Msa(pdbGenomeMsas);
		System.out.println("MSA:\t" + gp.getMsaId() + ":" + gp.getAaIdx());
	}
}
