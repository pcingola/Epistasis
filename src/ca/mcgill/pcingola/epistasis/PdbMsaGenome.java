package ca.mcgill.pcingola.epistasis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Stream;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

import ca.mcgill.mcb.pcingola.interval.Gene;
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
	public static final double PDB_RESOLUTION = 3.0; // PDB file resolution (in Angstrom)
	public static final String PDB_ORGANISM_COMMON = "HUMAN"; // PDB organism

	// Select ID
	public static final Function<IdMapperEntry, String> idme2id = ime -> ime.refSeqId;

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
		Timer.showStdErr("Checking PDB <-> Transcript sequences");

		// Initialize trancriptById
		if (trancriptById == null) {
			trancriptById = new HashMap<>();
			for (Gene g : config.getSnpEffectPredictor().getGenome().getGenes())
				for (Transcript tr : g) {
					String id = tr.getId();
					id = id.substring(0, id.indexOf('.'));
					trancriptById.put(id, tr);
				}
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
		String trIdsStr = IdMapper.ids(idMapper.getByPdbId(pdbId), idme2id);
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
					} else countMismatch++;
					sb.append(aa.getChemComp().getOne_letter_code());
				}

			// Only use mappings that have low error rate
			if (countMatch + countMismatch > 0) {
				double err = countMismatch / ((double) (countMatch + countMismatch));
				if (debug) System.out.println("\tChain: " + chain.getChainID() + "\terror: " + err + "\t" + sb);

				if (err < MAX_MISMATCH_RATE) {
					if (debug) Gpr.debug("Confirm transcript " + trId + "\terror: " + err);

					idmapsOri.stream() //
							.filter(idm -> trId.equals(idme2id.apply(idm)) && pdbId.equals(idm.pdbId)) //
							.findFirst() //
							.ifPresent(idm -> idmapsNew.add(idm));
				}
			}
		}

		if (verbose) idmapsNew.stream().forEach(i -> System.out.println("Confirmed IdMapping:\t" + i));

		return idmapsNew;
	}

	/**
	 * Distance analysis: Calculate amino acids closer than 'distanceThreshold' (within each protein) 
	 */
	public void distanceAnalysis(double distanceThreshold) {
		Timer.showStdErr("Distance analysis");
		// Get all results 
		PdbDistanceAnalysis pda = new PdbDistanceAnalysis(pdbDir, distanceThreshold, idMapper);
		List<DistanceResult> results = pda.run();

		// Map them to MSA
		Timer.showStdErr("Distance analysis: Find sequences from results");
		for (DistanceResult dres : results) {
			List<IdMapperEntry> idmes = idMapper.getByPdbId(dres.pdbId);

			// Find all transcripts
			for (IdMapperEntry idme : idmes) {
				String trid = idme2id.apply(idme);
				Transcript tr = trancriptById.get(trid);
				if (tr == null) throw new RuntimeException("Transcript '" + trid + "' not found. This should never happen!");

				// Find genomic position based on AA position
				int aa2pos[] = tr.aaNumber2Pos();
				if ((aa2pos.length <= dres.aaPos1) || (aa2pos.length <= dres.aaPos2)) {
					// System.out.println("\tPosition outside amino acid\tAA length: " + aa2pos.length + "\t" + dres);
					continue;
				}
				int pos1 = aa2pos[dres.aaPos1];
				int pos2 = aa2pos[dres.aaPos2];

				// Find sequences
				String seq1 = msas.findColumnSequence(trid, tr.getChromosomeName(), pos1);
				if (seq1 != null) System.out.println(dres + "\t" + tr.getChromosomeName() + ":" + pos1 + "\t" + seq1);
				String seq2 = msas.findColumnSequence(trid, tr.getChromosomeName(), pos2);
				if (seq2 != null) System.out.println(dres + "\t" + tr.getChromosomeName() + ":" + pos2 + "\t" + seq2);
			}
		}
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
