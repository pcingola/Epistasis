package ca.mcgill.pcingola.epistasis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.NextProt;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTreeAa;

/**
 * This class has information from
 * 		- PDB
 * 		- Genome (SnpEff)
 * 		- Multiple sequence alignment set
 * 		- IdMapper
 *
 * @author pcingola
 */
public class PdbGenomeMsas extends SnpEff {

	public static final double MAX_MISMATCH_RATE = 0.1;
	public static final double PDB_RESOLUTION = 3.0; // PDB file resolution (in Angstrom)
	public static final String PDB_ORGANISM_COMMON = "HUMAN"; // PDB organism
	public static int MAX_WARN = 20;

	int warn;
	public CountByType countMatch = new CountByType();
	CountByType countWarn = new CountByType();
	String genomeVer, pdbDir;
	HashMap<String, Transcript> trancriptById;
	IdMapper idMapper;
	LikelihoodTreeAa tree;
	MultipleSequenceAlignmentSet msas;
	PDBFileReader pdbreader;

	public PdbGenomeMsas(String configFile, String genomeVer, String pdbDir, MultipleSequenceAlignmentSet msas) {
		super(null);
		this.configFile = configFile;
		this.genomeVer = genomeVer;
		this.pdbDir = pdbDir;
		this.msas = msas;
	}

	/**
	 * Check that all protein sequences from genome match MSAs
	 */
	public void checkSequenceGenomeMsas() {
		trancriptById.keySet().stream().forEach(trid -> checkSequenceGenomeMsas(trid));
	}

	/**
	 * Check is the protein sequence from transcript 'trId' (genome) matches MSAs
	 * @return 'true' if protein sequences match
	 */
	boolean checkSequenceGenomeMsas(String trid) {
		Transcript tr = trancriptById.get(trid);
		if (tr == null) {
			warn("Transcript not found:", trid);
			return false;
		}

		// Get protein sequences
		String proteinTr = removeAaStop(tr.protein());
		String proteinMsa = removeAaStop(msas.findRowSequence(trid, tr.getChromosomeName()));

		if (!proteinTr.isEmpty() && proteinMsa != null) {
			boolean match = proteinTr.equals(proteinMsa);

			countMatch.inc("PROTEIN_MSA_VS_TR\t" + (match ? "OK" : "ERROR"));
			if (debug && !match) {
				System.err.println(trid + "\t" //
						+ (proteinTr.equals(proteinMsa) ? "OK" : "ERROR") //
						+ "\n\tPortein Tr  :\t" + proteinTr //
						+ "\n\tPortein MSA :\t" + proteinMsa //
						+ "\n");
			}

			return match;
		}

		return false;
	}

	/**
	 * Check that protein sequences form PDB matches sequences from Genome
	 * Return an IdMapped of confirmed entries (i.e. AA sequence matches between transcript and PDB)
	 */
	IdMapper checkSequencePdbGenome() {
		if (debug) System.err.println("Checking PDB <-> Transcript sequences\tdebug:" + debug);

		// Create a new IdMapper using only confirmed entries
		IdMapper idMapperConfirmed = new IdMapper();
		try {
			Files.list(Paths.get(pdbDir)) //
					.filter(s -> s.toString().endsWith(".pdb")) //
					.map(pf -> checkSequencePdbGenome(pf.toString())) //
					.forEach(ims -> idMapperConfirmed.addAll(ims));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		idMapper = idMapperConfirmed;
		return idMapperConfirmed;
	}

	/**
	 * Check that protein sequences form PDB (pdbFile) matches sequences from Genome
	 * Return a stream of maps that are confirmed (i.e. AA sequence matches between transcript and PDB)
	 */
	List<IdMapperEntry> checkSequencePdbGenome(String pdbFile) {
		Structure pdbStruct;
		try {
			pdbStruct = pdbreader.getStructure(pdbFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		String pdbId = pdbStruct.getPDBCode();

		// Get transcript IDs
		List<IdMapperEntry> idEntries = idMapper.getByPdbId(pdbId);
		String trIdsStr = IdMapper.ids(idEntries, IdMapperEntry.IDME_TO_REFSEQ);

		if (debug) {
			System.err.println(pdbId);
			System.err.println("\tEntries: ");
			if (idEntries != null) {
				idEntries.forEach(le -> System.err.println("\t\t" + le));
				System.err.println("\tTranscripts:\t" + trIdsStr);
			}
		}

		ArrayList<IdMapperEntry> list = new ArrayList<IdMapperEntry>();
		if (trIdsStr == null) return list;
		String trIds[] = trIdsStr.split(",");

		// Check idMaps. Only return those that match
		for (String trid : trIds)
			list.addAll(checkSequencePdbGenome(pdbId, pdbStruct, trid));

		return list;
	}

	/**
	 * Check that protein sequences match between PDB and Genome
	 * Return a list of maps that are confirmed (i.e. AA sequence matches between transcript and PDB)
	 * Note: Only part of the sequence usually matches
	 */
	List<IdMapperEntry> checkSequencePdbGenome(String pdbId, Structure pdbStruct, String trId) {
		if (debug) System.err.println("\nChecking " + trId + "\t<->\t" + pdbStruct.getPDBCode());
		List<IdMapperEntry> idmapsOri = idMapper.getByPdbId(pdbId);
		List<IdMapperEntry> idmapsNew = new ArrayList<>();

		// Transcript
		Transcript tr = trancriptById.get(trId);
		if (tr == null) return idmapsNew;
		String prot = tr.protein();
		if (debug) System.err.println("\tProtein: " + prot);

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
				if (debug) System.err.println("\tChain: " + chain.getChainID() + "\terror: " + err + "\t" + sb);

				if (err < MAX_MISMATCH_RATE) {
					if (debug) System.err.println("\t\tMapping OK    :\t" + trId + "\terror: " + err);

					int trAaLen = tr.protein().length();
					int pdbAaLen = chain.getAtomGroups("amino").size();

					idmapsOri.stream() //
							.filter(idm -> trId.equals(IdMapperEntry.IDME_TO_REFSEQ.apply(idm)) && pdbId.equals(idm.pdbId)) //
							.findFirst() //
							.ifPresent(i -> idmapsNew.add(i.cloneAndSet(chain.getChainID(), pdbAaLen, trAaLen)));
				} else if (debug) System.err.println("\t\tMapping ERROR :\t" + trId + "\terror: " + err);
			}
		}

		// Show all confirmed mappings
		idmapsNew.stream().forEach(i -> System.out.println(i));

		return idmapsNew;
	}

	/**
	 * Find a multiple sequence alignment based on a transcript ID and a genomic position
	 */
	public Triplet<String, String, Integer> findColumnSequence(String trid, int pos) {
		// Find all MSA
		List<MultipleSequenceAlignment> msaList = msas.getMsasByTrId(trid);
		if (msaList == null) return null;

		Transcript tr = trancriptById.get(trid);

		// Check all MSA
		for (MultipleSequenceAlignment msa : msaList) {
			// Different chromosome or position? Skip
			if (!msa.getChromo().equals(tr.getChromosomeName())) continue;
			if (pos < msa.getStart() || msa.getEnd() < pos) continue;

			// Find exon
			Exon exon = tr.findExon(pos);
			if (exon == null) {
				Gpr.debug("Cannot find exon for position " + pos + " in transcript " + tr.getId());
				return null;
			}

			// Find index
			int idxBase = tr.isStrandPlus() ? (pos - msa.getStart()) : (msa.getEnd() - pos);
			int idxAa = idxBase / 3;

			// WARNIGN: If exon frame is 1, the MSA has one additional AA (from the previous exon).
			//          I don't know why they do it this way...
			if (exon.getFrame() == 1) idxAa++;

			// Return column sequence
			return new Triplet<String, String, Integer>(msa.getColumnString(idxAa), msa.getId(), idxAa);
		}

		return null;
	}

	/**
	 * Load all data
	 */
	void initialize() {
		//---
		// Initialize SnpEff
		//---

		String argsSnpEff[] = { "eff", "-v", "-c", configFile, genomeVer };
		args = argsSnpEff;
		setGenomeVer(genomeVer);
		parseArgs(argsSnpEff);
		loadConfig();

		// Load SnpEff database
		if (genomeVer != null) loadDb();

		// Initialize trancriptById
		trancriptById = new HashMap<>();
		for (Gene g : config.getSnpEffectPredictor().getGenome().getGenes())
			for (Transcript tr : g) {
				String id = tr.getId();
				if (id.indexOf('.') > 0) id = id.substring(0, id.indexOf('.')); // When using RefSeq transcripts, we don't store sub-version number
				trancriptById.put(id, tr);
			}

		//---
		// Initialize reader
		//---
		pdbreader = new PDBFileReader();
	}

	/**
	 * Map 'DistanceResult' (Pdb coordinates) to MSA (genomic coordinates)
	 */
	public void mapToMsa(MultipleSequenceAlignmentSet msas, DistanceResult dres) {
		List<IdMapperEntry> idmes = idMapper.getByPdbId(dres.pdbId, dres.pdbChainId);
		if (idmes == null || idmes.isEmpty()) {
			warn("No mapping found for PdbId: ", "'" + dres.pdbId + "', chain '" + dres.pdbChainId + "'");
			return;
		}

		// Find all transcripts
		for (IdMapperEntry idme : idmes) {
			String trid = IdMapperEntry.IDME_TO_REFSEQ.apply(idme);

			// Transcript's protein doesn't match MSA's protein? Nothing to do
			if (!checkSequenceGenomeMsas(trid)) {
				countMatch.inc("_Total\tERROR\tTR-MSA");
				return;
			}

			Transcript tr = trancriptById.get(trid);
			if (tr == null) {
				warn("Transcript not found", trid);
				return;
			}

			// Find genomic position based on AA position
			int aa2pos[] = tr.aaNumber2Pos();
			if ((aa2pos.length <= dres.aaPos1) //
					|| (aa2pos.length <= dres.aaPos2) //
					|| (dres.aaPos1 < 0) //
					|| (dres.aaPos2 < 0) //
			) {
				// Position outside amino acid
				continue;
			}

			// Convert to genomic positions
			int pos1 = aa2pos[dres.aaPos1];
			int pos2 = aa2pos[dres.aaPos2];

			// Find sequences
			Triplet<String, String, Integer> res1 = findColumnSequence(trid, pos1);
			Triplet<String, String, Integer> res2 = findColumnSequence(trid, pos2);

			// Both sequences are available?
			if ((res1 != null) && (res2 != null)) {
				String seq1 = res1.a;
				String seq2 = res2.a;
				if ((seq1 != null) && (seq2 != null)) {

					Exon exon1 = tr.findExon(pos1);
					Exon exon2 = tr.findExon(pos2);

					// Count correct mappings
					boolean ok = (dres.aa1 == seq1.charAt(0)) && (dres.aa2 == seq2.charAt(0));
					String okStr = (ok ? "OK___" : "ERROR");
					String ok1Str = dres.aa1 != seq1.charAt(0) ? "ERROR" : "OK___";
					String ok2Str = dres.aa2 != seq2.charAt(0) ? "ERROR" : "OK___";
					countMatch.inc("_TOTAL_" + okStr);

					// Detailed counts for debugging
					if (debug) {
						countMatch.inc(dres.pdbId + "_" + ok1Str);
						countMatch.inc(dres.pdbId + "_" + ok2Str);
						countMatch.inc("_TOTAL_" + ok1Str + "_Strand:" + (tr.isStrandPlus() ? "+" : "-") + "_Frame:" + exon1.getFrame());
						countMatch.inc("_TOTAL_" + ok2Str + "_Strand:" + (tr.isStrandPlus() ? "+" : "-") + "_Frame:" + exon2.getFrame());
					}

					// Add information
					if (ok) {
						dres.chr1 = tr.getChromosomeName();
						dres.pos1 = pos1;
						dres.aaSeq1 = seq1;

						dres.chr2 = tr.getChromosomeName();
						dres.pos2 = pos2;
						dres.aaSeq2 = seq2;

						dres.transcriptId = tr.getId();

						dres.msa1 = res1.b;
						dres.msaIdx1 = res1.c;
						dres.msa2 = res2.b;
						dres.msaIdx2 = res2.c;
					} else {
						// Show mapping errors
						if (debug) System.err.println(ok1Str + " " + ok2Str //
								+ "\t" + dres.pdbId //
								+ "\t" + tr.getId() //
								+ "\t" + dres.distance //
								+ "\n\t" + dres.aa1 + "\t" + dres.aaPos1 + "\t" + tr.getChromosomeName() + ":" + pos1 + "\t" + exon1.getFrame() + "\t" + seq1 //
								+ "\n\t" + dres.aa2 + "\t" + dres.aaPos2 + "\t" + tr.getChromosomeName() + ":" + pos2 + "\t" + exon2.getFrame() + "\t" + seq2 //
						);
					}
				}
			}
		}
	}

	/**
	 * Create a marker encomapsing an amino acid (trId:aaIdx)
	 */
	protected Marker markerMsa(String trId, String chr, int start, int end, int aaIdx, char aaExpected) {
		// Find transcript and exon
		Transcript tr = trancriptById.get(trId);
		if (tr == null) {
			Gpr.debug("TR NOT FOUND: " + trId);
			return null;
		}
		Exon ex = tr.findExon(start);
		if (ex == null) return null;

		// Calculate start position
		int startPos;
		int fr = 0;
		if (ex.getFrame() != 0) {
			if (ex.isStrandPlus()) {
				aaIdx--;
				if (ex.getFrame() == 2) aaIdx++; // I don't know why UCSC numbers the AA differentlt when frame is 2
				fr = 3 - ex.getFrame(); // Offset based on frame
			} else {
				aaIdx--;
				if (ex.getFrame() == 2) aaIdx++; // I don't know why UCSC numbers the AA differentlt when frame is 2
				fr = 3 - ex.getFrame(); // Offset based on frame
			}
		}

		// Find AA start position
		if (ex.isStrandPlus()) {
			int exStart = Math.max(start, tr.getCdsStart());
			startPos = exStart + (aaIdx * 3 + fr);
		} else {
			int exEnd = Math.min(end, tr.getCdsStart());
			startPos = exEnd - (aaIdx * 3 + fr);
		}

		// Get position within CDS
		int cdsBase = tr.baseNumberCds(startPos, false);
		int cds2pos[] = tr.baseNumberCds2Pos();
		if ((ex.isStrandPlus() && (startPos < ex.getStart())) //
				|| (ex.isStrandMinus() && (startPos > ex.getEnd()))) {
			// If the position is outside the exon, then we must jump to previous exon
			startPos = cds2pos[cdsBase - ex.getFrame()];
			cdsBase = tr.baseNumberCds(startPos, true);
		}

		//---
		// Sanity check: Make sure that AA matches between transcript model and MSA data from 'genes likelihood' file
		//---
		String id = chr + ":" + start + "-" + end + "[" + aaExpected + "]";

		// Extract codon
		String cdsSeq = tr.cds();
		String codonStr = cdsSeq.substring(cdsBase, cdsBase + 3);
		String aa = genome.codonTable().aa(codonStr);

		if (aa.equals("" + aaExpected)) {

			if (debug) Gpr.debug("OK: " + id + " : " + aa);
		} else {
			if (debug) Gpr.debug("Entry ID     : " + id //
					+ "\ntr ID        : " + trId + ", chr: " + chr + ", start: " + start + ", end: " + end + ", idx: " + aaIdx + ", fr: " + fr//
					+ "\nTranscript : " + tr //
					+ "\nExon       : " + ex //
					+ "\nStart pos: " + startPos //
					+ "\nCodon    : " + codonStr + ", aa (real): " + aa + ", aa (exp): " + aaExpected //
			);
			return null;
		}

		//---
		// Create marker
		// Important: The marker has ALL bases in the codon.
		//            For instance, is the codon is split between two exons, the
		//            marker will contain the intron
		//---
		int markerStart, markerEnd;
		if (tr.isStrandPlus()) {
			markerStart = cds2pos[cdsBase];
			markerEnd = cds2pos[cdsBase + 2];
		} else {
			markerStart = cds2pos[cdsBase + 2];
			markerEnd = cds2pos[cdsBase];
		}

		Chromosome chromo = genome.getChromosome(chr);
		Marker marker = new Marker(chromo, markerStart, markerEnd, ex.isStrandMinus(), id);
		return marker;
	}

	/**
	 * Find nextprot annotations
	 */
	public void nextProt(DistanceResult d) {
		d.annotations1 = nextProt(d.chr1, d.pos1, d.transcriptId);
		d.annotations2 = nextProt(d.chr2, d.pos2, d.transcriptId);
		System.out.println(d);
	}

	/**
	 * Find nextprot annotations
	 */
	public String nextProt(String chr, int pos, String refSeqId) {
		if (refSeqId.indexOf('.') > 0) refSeqId = refSeqId.substring(0, refSeqId.indexOf('.'));

		// Create a set of transcript IDs
		List<IdMapperEntry> idEntries = idMapper.getByRefSeqId(refSeqId);
		String trIdsStr = IdMapper.ids(idEntries, IdMapperEntry.IDME_TO_ENSEMBLID);
		if (trIdsStr == null) return "";

		HashSet<String> trIds = new HashSet<String>();
		Arrays.stream(trIdsStr.split(",")).forEach(id -> trIds.add(id));

		// Find all NextProt entries matching any transcript ID
		Marker m = new Marker(config.getGenome().getChromosome(chr), pos, pos, false, "");
		Markers results = config.getSnpEffectPredictor().query(m);

		return results.stream() //
				.filter(r -> r instanceof NextProt && trIds.contains(((Transcript) r.getParent().getParent()).getId())) //
				.map(r -> r.getId()) //
				.sorted() //
				.distinct() //
				.collect(Collectors.joining(";") //
				);
	}

	/**
	 * Remove 'stop' AA form sequence
	 */
	String removeAaStop(String seq) {
		if (seq == null) return null;
		if (seq.endsWith("-") || seq.endsWith("*")) return seq.substring(0, seq.length() - 1);
		return seq;
	}

	public void resetStats() {
		countMatch = new CountByType();
	}

	public void setIdMapper(IdMapper idMapper) {
		this.idMapper = idMapper;
	}

	public void setNextProt(boolean nextProt) {
		this.nextProt = nextProt;
	}

	public void setTree(LikelihoodTreeAa tree) {
		this.tree = tree;
	}

	void warn(String warningType, String warning) {
		countWarn.inc(warningType);
		if (countWarn.getCount(warningType) < MAX_WARN) System.err.println("WARNING: " + warningType + " " + warning);
	}

}
