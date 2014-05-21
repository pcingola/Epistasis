package ca.mcgill.pcingola.epistasis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTree;

public class Zzz extends SnpEff {

	public static final double MAX_MISMATCH_RATE = 0.1;

	public static void main(String[] args) {
		// SnpEff database
		String genome = "testHg3771Chr1";
		// String argsSnpEff[] = { "eff", "-v", "-nextProt", "-c", Gpr.HOME + "/snpEff/" + Config.DEFAULT_CONFIG_FILE, genome };
		String argsSnpEff[] = { "eff", "-v", "-c", Gpr.HOME + "/snpEff/" + Config.DEFAULT_CONFIG_FILE, genome };

		Zzz zzz = new Zzz(argsSnpEff);
		zzz.setGenomeVer(genome);
		zzz.parseArgs(argsSnpEff);
		zzz.loadConfig();
		zzz.loadAll();
		zzz.loadDb();

		zzz.checkCoordinates();
	}

	String pdbId = "1AN4";
	String pdbFile = Gpr.HOME + "/snpEff/db/pdb/" + pdbId.toLowerCase() + ".pdb";
	String phyloFile = Gpr.HOME + "/snpEff/db/multiz100way/hg19.100way.nh";
	String multAlignFile = Gpr.HOME + "/snpEff/db/multiz100way/refGene.exonAA.head.fa";
	String idMapFile = Gpr.HOME + "/snpEff/db/multiz100way/idMap_ensemblId_refseq_pdbId.txt";
	IdMapper idMapper;
	LikelihoodTree tree;
	MultipleSequenceAlignmentSet msas;
	Structure structure;
	HashMap<String, Transcript> trancriptById;

	public Zzz(String[] args) {
		super(args);
	}

	/**
	 * Load files
	 */
	void loadAll() {
		// Pdb file
		Timer.showStdErr("Loading pdb file " + pdbFile);
		PDBFileReader pdbreader = new PDBFileReader();
		try {
			structure = pdbreader.getStructure(pdbFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		System.out.println(structure);

		// Id Map
		Timer.showStdErr("Loading id maps " + idMapFile);
		idMapper = new IdMapper(idMapFile);

		// Load: tree
		Timer.showStdErr("Loading phylogenetic tree from " + phyloFile);
		tree = new LikelihoodTree();
		tree.load(phyloFile);
		int numAligns = tree.childNames().size();

		// Load: MSA
		Timer.showStdErr("Loading " + numAligns + " way multiple alignment from " + multAlignFile);
		msas = new MultipleSequenceAlignmentSet(multAlignFile, numAligns);
		msas.load();
	}

	/**
	 * Map coordinates from Pdb to Transcript
	 */
	void checkCoordinates() {
		Timer.showStdErr("Mapping coordinates");

		// Initialize trancriptById
		if (trancriptById == null) {
			trancriptById = new HashMap<>();
			config.getSnpEffectPredictor().getGenome().getGenes().forEach( //
					g -> g.subintervals().forEach( //
							t -> trancriptById.put(t.getId(), t)) //
					);
		}

		// Get trancsript IDs
		String trIdsStr = IdMapper.trIds(idMapper.getByPdbId(pdbId));
		System.out.println("Mapping Pdb ID '" + pdbId + "' to transcripts: " + trIdsStr);
		String trIds[] = trIdsStr.split(",");

		// Mapping coordinates
		HashSet<IdMapperEntry> confirmed = new HashSet<>();
		Arrays.stream(trIds).forEach(id -> confirmed.addAll(checkCoordinates(pdbId, structure, id)));
		Gpr.debug("Confirmed: ");
		confirmed.stream().forEach(System.out::println);
	}

	/**
	 * Return a list of maps that are confirmed (i.e. AA sequence matches between transcript and PDB)
	 * Note: Only part of the sequence usually matches
	 */
	List<IdMapperEntry> checkCoordinates(String pdbId, Structure pdbStruct, String trId) {
		if (debug) System.out.println("Checking " + trId + "\t<->\t" + pdbStruct.getPDBCode());
		List<IdMapperEntry> idmapsOri = idMapper.getByPdbId(pdbId);
		List<IdMapperEntry> idmapsNew = new ArrayList<>();

		// Transcript
		Transcript tr = trancriptById.get(trId);
		if (tr == null) return idmapsNew;
		String prot = tr.protein();
		if (debug) System.out.println("\tProtein: " + prot);

		// Compare to PDB structure
		for (Chain chain : structure.getChains()) {
			// Compare sequence to each AA-Chain
			StringBuilder sb = new StringBuilder();
			int countMatch = 0, countMismatch = 0;

			// Count differences
			for (Group group : chain.getAtomGroups())
				if (group instanceof AminoAcid) {
					AminoAcid aa = (AminoAcid) group;
					int aaPos = aa.getResidueNumber().getSeqNum() - 1;
					char aaLetter = aa.getChemComp().getOne_letter_code().charAt(0);
					if (prot.length() > aaPos) {
						char trAaLetter = prot.charAt(aaPos);

						if (aaLetter == trAaLetter) countMatch++;
						else countMismatch++;

						if (debug) System.out.println("\t\t" + aaPos + "\t" + aaLetter + "\t" + trAaLetter + "\t" + (aaLetter != trAaLetter ? "*" : ""));
					} else countMismatch++;
					sb.append(aa.getChemComp().getOne_letter_code());
				}

			// Only use mappings that have low error rate
			if (countMatch + countMismatch > 0) {
				double err = countMismatch / ((double) (countMatch + countMismatch));
				if (debug) System.out.println("\tChain: " + chain.getChainID() + "\terror: " + err + "\t" + sb);

				if (err < MAX_MISMATCH_RATE) {
					if (debug) Gpr.debug("Confirm transcript " + trId);

					idmapsOri.stream() //
							.filter(idm -> trId.equals(idm.trId) && pdbId.equals(idm.pdbId)) //
							.findFirst() //
							.ifPresent(idm -> idmapsNew.add(idm));
				}
			}
		}

		return idmapsNew;
	}
}
