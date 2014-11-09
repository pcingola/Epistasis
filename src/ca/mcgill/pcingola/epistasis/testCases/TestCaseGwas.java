package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;

import org.junit.Assert;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.GenotypePos;
import ca.mcgill.pcingola.epistasis.gwas.GwasEpistasis;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.pdb.PdbGenomeMsas;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTreeAa;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseGwas extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	public void test_01_Gwas_Map() {
		Gpr.debug("Test");

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String genesLikeFile = "test/NM_001438.txt";
		String vcfFile = ""; // It doesn't matter, it is not used
		String phenoFile = ""; // It doesn't matter, it is not used

		GwasEpistasis gwasEpistasis = new GwasEpistasis(configFile, genome, genesLikeFile, vcfFile, phenoFile);
		gwasEpistasis.setDebug(debug);
		gwasEpistasis.initialize();
		gwasEpistasis.readGenesLogLikelihood();

		// Check that all mappings are OK
		Assert.assertEquals(0, gwasEpistasis.getCountErr());
		Assert.assertEquals(334, gwasEpistasis.getCountOk());
	}

	public void test_02_Gwas_Map() {
		Gpr.debug("Test");

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String genesLikeFile = "test/NM_021969.txt";
		String vcfFile = ""; // It doesn't matter, it is not used
		String phenoFile = ""; // It doesn't matter, it is not used

		GwasEpistasis gwasEpistasis = new GwasEpistasis(configFile, genome, genesLikeFile, vcfFile, phenoFile);
		gwasEpistasis.setDebug(debug);
		gwasEpistasis.initialize();
		gwasEpistasis.readGenesLogLikelihood();

		// Check that all mappings are OK
		Assert.assertEquals(0, gwasEpistasis.getCountErr());
		Assert.assertEquals(458, gwasEpistasis.getCountOk());
	}

	public void test_03_Gwas_Map() {
		Gpr.debug("Test");

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String genesLikeFile = "test/NM_004905.txt";
		String vcfFile = ""; // It doesn't matter, it is not used
		String phenoFile = ""; // It doesn't matter, it is not used

		GwasEpistasis gwasEpistasis = new GwasEpistasis(configFile, genome, genesLikeFile, vcfFile, phenoFile);
		gwasEpistasis.setDebug(debug);
		gwasEpistasis.initialize();
		gwasEpistasis.readGenesLogLikelihood();

		// Check that all mappings are OK
		Assert.assertEquals(0, gwasEpistasis.getCountErr());
		Assert.assertEquals(298, gwasEpistasis.getCountOk());
	}

	public void test_04_Gwas_Map() {
		Gpr.debug("Test");

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String genesLikeFile = "test/gwas_map_test_chr1.txt";
		String vcfFile = ""; // It doesn't matter, it is not used
		String phenoFile = ""; // It doesn't matter, it is not used

		GwasEpistasis gwasEpistasis = new GwasEpistasis(configFile, genome, genesLikeFile, vcfFile, phenoFile);
		gwasEpistasis.setDebug(debug);
		gwasEpistasis.initialize();
		gwasEpistasis.readGenesLogLikelihood();

		// Check that all mappings are OK
		Assert.assertEquals(0, gwasEpistasis.getCountErr());
		Assert.assertEquals(13280, gwasEpistasis.getCountOk());
	}

	public void test_05_Gwas_Map_RoundTrip() {
		// Create a atest to map using
		//	i) Select <msaId, aaIdx>
		//	ii) Map it to genomic coordinate
		// 	iii) Map genomic coordinates back to msaId:aaIdx
		//  iv) Check that <mdsId, aaIdx> are recovered correctly

		String configFile = Gpr.HOME + "/snpEff/snpEff.config";
		String genome = "testHg19Chr1";
		String phyloFileName = "data/hg19.100way.nh";
		String msasFile = "data/msa_test.fa.gz";
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
		int n = 0;
		for (MultipleSequenceAlignment msa : msas) {
			String trId = msa.getTranscriptId();

			// Sanity check: Do protein sequences match? (transcript vs MSA)
			if (!pdbGenomeMsas.checkSequenceGenomeMsas(trId)) {
				if (pdbGenomeMsas.getTranscript(trId) != null) {
					if (debug) Gpr.debug("AA Sequences differ:" //
							+ "\n\tMSA        : " + msas.rowSequence(trId) //
							+ "\n\tTranscript : " + pdbGenomeMsas.getTranscript(trId).protein() //
					);
				}
				continue;
			}

			// Step i: Select all <msaId, aaIdx>
			for (int aaIdx = 0; aaIdx < msa.getAaSeqLen(); aaIdx++) {

				// Step ii: Map it to genomic coordinate
				GenotypePos gp = new GenotypePos(msa.getId(), aaIdx);
				String err = gp.mapMsa2GenomicErr(pdbGenomeMsas);
				if (err == null) {
					if (verbose) Gpr.debug(n + "\tOK\t" + gp.getMsaId() + ":" + gp.getAaIdx() + "\t" + gp.getChromosomeName() + ":" + gp.getStart() + "-" + gp.getEnd());
				} else {
					if (verbose) Gpr.debug(n + "\tNO\t" + gp.getMsaId() + ":" + gp.getAaIdx() + "\t" + err);

					String acceptedErr = "Transcript '" + trId + "' not found";
					if (!err.equals(acceptedErr)) throw new RuntimeException("Unacceptable error condition: " + err);

					continue;
				}

				// Step iii:  Map genomic coordinates back to msaId:aaIdx
				GenotypePos gp2 = new GenotypePos(gp.getParent(), gp.getStart(), gp.getId());
				if (!gp2.mapGenomic2Msa(pdbGenomeMsas, trId)) continue;

				// Step iv: Check that coordinates are mapped back correctly
				String msaIdxOri = gp.getMsaId() + ":" + gp.getAaIdx();
				String msaIdxRecover = gp2.getMsaId() + ":" + gp2.getAaIdx();

				if (!msaIdxOri.equals(msaIdxRecover) || verbose) Gpr.debug(n //
						+ "\t" + msaIdxOri //
						+ "\t" + gp.getChromosomeName() + ":" + gp.getStart() + "-" + gp.getEnd() //
						+ "\t" + msaIdxRecover //
				);

				Assert.assertEquals(msaIdxOri, msaIdxRecover);

				n++;
				if (!verbose) Gpr.showMark(n, 100);
			}
		}
	}
}
