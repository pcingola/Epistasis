package ca.mcgill.pcingola.epistasis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Stream;

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

/**
 * This class has information from
 * 		- PDB
 * 		- Multiple Sequence alignment
 * 		- Genome (SnpEff)
 * 		- Phylogeny tree
 * 		- IdMapper
 *
 * @author pcingola
 *
 */
public class PdbMsaGenome extends SnpEff {

	public static final double MAX_MISMATCH_RATE = 0.1;

	public static void main(String[] args) {
		String genome = "testHg3771Chr1";
		PdbMsaGenome zzz = new PdbMsaGenome( //
				Gpr.HOME + "/snpEff/" + Config.DEFAULT_CONFIG_FILE //
				, genome //
				, Gpr.HOME + "/snpEff/db/pdb" //
				, Gpr.HOME + "/snpEff/db/multiz100way/hg19.100way.nh" //
				, Gpr.HOME + "/snpEff/db/multiz100way/refGene.exonAA.head.fa" //
				, Gpr.HOME + "/snpEff/db/multiz100way/idMap_ensemblId_refseq_pdbId.txt" //
		);

		zzz.initialize();
		zzz.checkCoordinates();
	}

	public static final double PDB_RESOLUTION = 3.0; // PDB file resolution (in Angstrom)
	public static final String PDB_ORGANISM_COMMON = "HUMAN"; // PDB organism

	String genome, pdbDir, phyloFile, multAlignFile, idMapFile;
	IdMapper idMapper;
	LikelihoodTree tree;
	MultipleSequenceAlignmentSet msas;
	HashMap<String, Transcript> trancriptById;
	PDBFileReader pdbreader;

	public PdbMsaGenome(String args[]) {
		this(args[0], args[1], args[2], args[3], args[4], args[5]);
	}

	public PdbMsaGenome(String configFile, String genome, String pdbDir, String phyloFile, String multAlignFile, String idMapFile) {
		super(null);
		this.configFile = configFile;
		this.genome = genome;
		this.pdbDir = pdbDir;
		this.phyloFile = phyloFile;
		this.multAlignFile = multAlignFile;
		this.idMapFile = idMapFile;
	}

	/**
	 * Check mapping.
	 * Return an IdMapped of confirmed entries (i.e. AA sequence matches between transcript and PDB)
	 */
	IdMapper checkCoordinates() {
		Timer.showStdErr("Mapping coordinates");

		// Initialize trancriptById
		if (trancriptById == null) {
			trancriptById = new HashMap<>();
			config.getSnpEffectPredictor().getGenome().getGenes().forEach( //
					g -> g.subintervals().forEach( //
							t -> trancriptById.put(t.getId(), t)) //
					);
		}

		// Create a new IdMapper using only confirmed entries
		IdMapper idMapperConfirmed = new IdMapper();
		try {
			Files.list(Paths.get(pdbDir)) //
					.filter(s -> s.toString().endsWith(".pdb")) //
					.map(pf -> checkCoordinates(pf.toString())) //
					.flatMap(s -> s) //
					.forEach(im -> idMapperConfirmed.add(im));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		idMapper = idMapperConfirmed;
		return idMapperConfirmed;
	}

	/**
	 * Check mapping.
	 * Return a stream of maps that are confirmed (i.e. AA sequence matches between transcript and PDB)
	 */
	Stream<IdMapperEntry> checkCoordinates(String pdbFile) {
		Structure pdbStruct;
		try {
			pdbStruct = pdbreader.getStructure(pdbFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		String pdbId = pdbStruct.getPDBCode();

		// Get trancsript IDs
		String trIdsStr = IdMapper.trIds(idMapper.getByPdbId(pdbId));
		String trIds[] = trIdsStr.split(",");

		// Check idMaps. Only return those that match
		return Arrays.stream(trIds) //
				.map(id -> checkCoordinates(pdbId, pdbStruct, id)) //
				.flatMap(l -> l.stream());
	}

	/**
	 * Check mapping.
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

		// Filter PDB structure
		// Within resolution limits? => Process
		double res = pdbStruct.getPDBHeader().getResolution();
		if (res > PDB_RESOLUTION) return idmapsNew;

		// Compare to PDB structure
		for (Chain chain : pdbStruct.getChains()) {
			// Compare sequence to each AA-Chain
			StringBuilder sb = new StringBuilder();
			int countMatch = 0, countMismatch = 0;

			// Count differences
			for (Group group : chain.getAtomGroups())
				if (group instanceof AminoAcid) {
					AminoAcid aa = (AminoAcid) group;
					int aaPos = aa.getResidueNumber().getSeqNum() - 1;
					if (aaPos < 0) continue; // I don't know why some PDB coordinates are negative...

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

		if (verbose) {
			System.out.println("Mapping Pdb ID\t" + pdbId);
			idmapsNew.stream().forEach(i -> System.out.println("Confirmed IdMapping:\t" + i));
		}

		return idmapsNew;
	}

	/**
	 * Load all data
	 */
	void initialize() {
		// Set up SnpEff arguments
		String argsSnpEff[] = { "eff", "-v", "-c", configFile, genome };
		args = argsSnpEff;
		setGenomeVer(genome);
		parseArgs(argsSnpEff);
		loadConfig();

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

		// Initialize reader
		pdbreader = new PDBFileReader();

		// Load SnpEff database
		loadDb();
	}
}
