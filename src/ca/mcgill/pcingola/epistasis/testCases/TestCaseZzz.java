package ca.mcgill.pcingola.epistasis.testCases;

import java.util.Random;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.GenotypePos;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignment;
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

	public void test_05_Gwas_Map_RoundTrip() {
		// Create a atest to map using
		//	i) Select a random <msaId, aaIdx>
		//	ii) Map it to genomic coordinate 
		// 	iii) Map genomic coordinates back to msaId:aaIdx
		//  iv) Check that <mdsId, aaIdx> are recovered correctly

		//---
		// Initialize and load data
		//---
		Random random = new Random(20141006);

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String phyloFileName = "data/hg19.100way.nh";
		// String msasFile = "data/msa_test.fa.gz";
		String msasFile = "data/hm1.fa";
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

		//---
		// Test N times
		//---
		int N = 1;

		for (int i = 0; i < N;) {
			// Step i: Select a random <msaId, aaIdx>
			MultipleSequenceAlignment msa = msas.rand(random);
			int aaIdx = random.nextInt(msa.getAaSeqLen());
			String trId = msa.getTranscriptId();

			// Step ii: Map it to genomic coordinate
			GenotypePos gp = new GenotypePos(msa.getId(), aaIdx);
			String err = gp.mapMsa2GenomicErr(pdbGenomeMsas);
			if (err == null) {
				if (verbose) Gpr.debug(i + "\tOK\t" + gp.getMsaId() + ":" + gp.getAaIdx() + "\t" + gp.getChromosomeName() + ":" + gp.getStart() + "-" + gp.getEnd());
			} else {
				if (verbose) Gpr.debug(i + "\tNO\t" + gp.getMsaId() + ":" + gp.getAaIdx() + "\t" + err);
				continue;
			}

			// Step iii:  Map genomic coordinates back to msaId:aaIdx
			GenotypePos gp2 = new GenotypePos(gp.getParent(), gp.getStart(), gp.getId());
			if (!gp2.mapGenomic2Msa(pdbGenomeMsas, trId)) continue;

			// Step iv: Check that coordinates are mapped back correctly
			String msaIdxOri = gp.getMsaId() + ":" + gp.getAaIdx();
			String msaIdxRecover = gp2.getMsaId() + ":" + gp2.getAaIdx();
			Gpr.debug(i + "\t" + (msaIdxOri.equals(msaIdxRecover)) //
					+ "\t" + msaIdxOri //
					+ "\t" + gp.getChromosomeName() + ":" + gp.getStart() + "-" + gp.getEnd() //
					+ "\t" + msaIdxRecover //
			);

			Assert.assertEquals(msaIdxOri, msaIdxRecover);
			i++;
		}

		throw new RuntimeException("Unimplemented test!");
	}
	//	public void test_06_Map_InDels() {
	//		throw new RuntimeException("Unimplemented test!");
	//	}

}
