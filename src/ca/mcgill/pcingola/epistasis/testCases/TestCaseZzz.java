package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
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

	public GenotypePos mapToMsa(String genome, String msasFile, String chr, int pos) {
		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String phyloFileName = "data/hg19.100way.nh";
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

		Chromosome chromo = pdbGenomeMsas.getConfig().getGenome().getOrCreateChromosome(chr);
		GenotypePos gp = new GenotypePos(chromo, pos - 1, chr + ":" + pos);
		gp.mapGenomic2Msa(pdbGenomeMsas);

		return gp;
	}

	public void test_zzz() {
		String genome = "testHg19Chr1";
		String msasFile = "test/NM_178229.fa";
		GenotypePos gp = mapToMsa(genome, msasFile, "1", 156526325);
		Assert.assertTrue(gp.getMsaId() != null);
		Timer.showStdErr("MSA:\t" + gp.getMsaId() + ":" + gp.getAaIdx());
	}
}
