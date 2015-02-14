package ca.mcgill.pcingola.epistasis.pdb;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.NextProt;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.epistasis.IdMapper;
import ca.mcgill.pcingola.epistasis.IdMapperEntry;
import ca.mcgill.pcingola.epistasis.coordinates.GenomicCoordinates;
import ca.mcgill.pcingola.epistasis.coordinates.MsaCoordinates;
import ca.mcgill.pcingola.epistasis.coordinates.PdbCoordinate;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;
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
	Map<String, Transcript> trancriptById;
	Set<String> trancriptsChecked;
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
		trancriptsChecked = new HashSet<String>();
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
	public boolean checkSequenceGenomeMsas(String trid) {
		Transcript tr = trancriptById.get(trid);
		if (tr == null) {
			warn("Transcript not found:", trid);
			return false;
		}

		// Already checked?
		if (trancriptsChecked.contains(trid)) return true;

		// Get protein sequences
		String proteinTr = removeAaStop(tr.protein());
		String proteinMsa = removeAaStop(msas.rowSequence(trid, tr.getChromosomeName()));

		if (!proteinTr.isEmpty() && proteinMsa != null) {
			boolean match = proteinTr.equals(proteinMsa);
			if (match) trancriptsChecked.add(trid);

			if (debug && !match) {
				countMatch.inc("PROTEIN_MSA_VS_TR\t" + (match ? "OK" : "ERROR"));
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
	public IdMapper checkSequencePdbGenome() {
		if (debug) System.err.println("Checking PDB <-> Transcript sequences\tdebug:" + debug);

		// Create a new IdMapper using only confirmed entries
		IdMapper idMapperConfirmed = new IdMapper();
		try {
			Files.list(Paths.get(pdbDir)) //
					.filter(s -> s.toString().endsWith(".pdb")) //
					.map(pf -> checkSequencePdbGenome(pf.toString(), true)) //
					.forEach(ims -> idMapperConfirmed.addAll(ims));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		idMapper = idMapperConfirmed;
		return idMapperConfirmed;
	}

	/**
	 * Check that protein sequences form PDB (pdbFile) matches sequences from Genome
	 * Return a list of maps that are confirmed (i.e. AA sequence matches between transcript and PDB)
	 */
	public List<IdMapperEntry> checkSequencePdbGenome(String pdbFile, boolean show) {
		Structure pdbStruct = readPdbFile(pdbFile);
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
			list.addAll(checkSequencePdbGenome(pdbId, pdbStruct, trid, show));

		return list;
	}

	/**
	 * Check that protein sequences match between PDB and Genome
	 * Return a list of maps that are confirmed (i.e. AA sequence matches between transcript and PDB)
	 * Note: Only part of the sequence usually matches
	 */
	List<IdMapperEntry> checkSequencePdbGenome(String pdbId, Structure pdbStruct, String trId, boolean show) {
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
		if (show) idmapsNew.stream().forEach(i -> System.out.println(i));

		return idmapsNew;
	}

	public MultipleSequenceAlignmentSet getMsas() {
		return msas;
	}

	public String getPdbDir() {
		return pdbDir;
	}

	public Transcript getTranscript(String trId) {
		return trancriptById.get(trId);
	}

	/**
	 * Load all data
	 */
	public void initialize() {
		//---
		// Initialize SnpEff
		//---
		String argsSnpEff[] = { "eff", "-c", configFile, genomeVer };
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

				if (trancriptById.containsKey(id)) {
					// There is already a transcript?
					// Favor shorter chromosome names. E.g.: 'chr6' is better than 'chr6_cox_hap2'
					String chrPrev = trancriptById.get(id).getChromosomeName();
					String chr = tr.getChromosomeName();

					if (chr.length() < chrPrev.length()) trancriptById.put(id, tr);
				} else {
					// Transcript not present: Add it
					trancriptById.put(id, tr);
				}
			}

		//---
		// Initialize reader
		//---
		pdbreader = new PDBFileReader();
	}

	/**
	 * Map a genomic position to an MSA index (given the Transcript and the MSA)
	 * Most of the times only aaIdx within 'msa' is calculated. On some
	 * rare border conditions", 'msa' changes (e.g. when 'pos' maps to
	 * the last AA and the next exon has frame=1). This case is a quirk
	 * on how UCSC numbers AA in their multiple sequence alignments.

	 * @return false on failure
	 */
	public MsaCoordinates mapMsaTrPos2AaIdx(MultipleSequenceAlignment msa, Transcript tr, int pos) {
		if (tr == null) return null;

		// Check all MSA
		// Different chromosome or position? Skip
		if (!msa.intersects(tr)) return null;

		// Find exon
		Exon exon = tr.findExon(pos);
		if (exon == null) {
			Gpr.debug("Cannot find exon for position " + pos + " in transcript " + tr.getId());
			return null;
		}

		// Find index
		int idxBase = tr.isStrandPlus() ? (pos - msa.getStart()) : (msa.getEnd() - pos);
		int idxAa = idxBase / 3;

		// WARNIGN: If exon frame is 2, the MSA has one additional AA (from the previous exon).
		//          I don't know why they do it this way...
		if (exon.getFrame() == 2) {
			if (idxBase < 1) idxAa = 0; // First two bases are AA number zero
			else idxAa++; // Other bases are AA number 1 and on
		}

		// Out of range
		if (idxAa >= msa.getAaSeqLen()) {
			// This can happen when a base maps to the LAST amino acid in an exon.
			// If the next exons has 'frame=2', then that last AA is 'pushed' to the
			// next exon (I don't know why UCSC does this complicated mapping between
			// AA and bases in their MSAs). So, we have to move the mapping to
			// the first AA in the next exon.
			MultipleSequenceAlignment msaNext = getMsas().findNextExon(msa); // Find MSA for the exon following 'msa'
			if (msaNext == null) {
				Gpr.debug("ERROR: Cannot find 'next' exon for MSA '" + msa.getId() + "'");
				return null; // Cannot find next exon (something went wrong)
			}
			msa = msaNext;
			idxAa = 0; // First amino acid
		}

		// We are done: Set parameters
		return new MsaCoordinates(msa.getId(), idxAa);
	}

	/**
	 * Map 'DistanceResult' (Pdb coordinates) to MSA (genomic coordinates)
	 */
	public void mapToMsa(DistanceResult dres) {
		// Find trancript IDs using PDB ids
		List<IdMapperEntry> idmes = idMapper.getByPdbId(dres.pdbId, dres.pdbChainId);
		if (idmes == null || idmes.isEmpty()) {
			warn("No mapping found for PdbId: ", "'" + dres.pdbId + "', chain '" + dres.pdbChainId + "'");
			return;
		}

		// Find all transcripts, then map <tr, pos> to <msaId, aaIdx>
		for (IdMapperEntry idme : idmes) {
			//---
			// Map <pdbId, aaPos> to <transcript, pos>
			//---
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

			//---
			// Map <transcript, pos> to <msaId, aaIdx>
			//---

			// Find MSA and aaIdx
			GenomicCoordinates gp1 = new GenomicCoordinates(), gp2 = new GenomicCoordinates();
			gp1.mapTrPos2MsaIdx(this, trid, pos1);
			gp2.mapTrPos2MsaIdx(this, trid, pos2);

			// Both mappings are available?
			if (gp1.hasMsaInfo() && gp2.hasMsaInfo()) {

				// Sanity check: Find sequences and check against PDB data
				String seq1 = msas.getMsa(gp1.getMsaId()).getColumnString(gp1.getAaIdx());
				String seq2 = msas.getMsa(gp2.getMsaId()).getColumnString(gp2.getAaIdx());

				if ((seq1 != null) && (seq2 != null)) {
					Exon exon1 = null, exon2 = null;
					// Does transcript AA sequence match PDB Aa sequence?
					boolean ok = (dres.aa1 == seq1.charAt(0)) && (dres.aa2 == seq2.charAt(0));

					// Count correct mappings
					String okStr = (ok ? "OK___" : "ERROR");
					String ok1Str = dres.aa1 != seq1.charAt(0) ? "ERROR" : "OK___";
					String ok2Str = dres.aa2 != seq2.charAt(0) ? "ERROR" : "OK___";
					countMatch.inc("_TOTAL_" + okStr);

					// Detailed counts for debugging
					if (debug) {
						exon1 = tr.findExon(pos1);
						exon2 = tr.findExon(pos2);
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

						dres.msa1 = gp1.getMsaId();
						dres.msaIdx1 = gp1.getAaIdx();
						dres.msa2 = gp2.getMsaId();
						dres.msaIdx2 = gp2.getAaIdx();
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
	 * Map 'DistanceResult' (Pdb coordinates) to MSA (genomic coordinates)
	 */
	public MsaCoordinates mapToMsa(PdbCoordinate pdb) {
		// Find trancript IDs using PDB ids
		List<IdMapperEntry> idmes = idMapper.getByPdbId(pdb.pdbId, pdb.pdbChainId);
		if (idmes == null || idmes.isEmpty()) {
			warn("No mapping found for PdbId: ", "'" + pdb.pdbId + "', chain '" + pdb.pdbChainId + "'");
			return null;
		}

		// Find all transcripts, then map <tr, pos> to <msaId, aaIdx>
		for (IdMapperEntry idme : idmes) {
			//---
			// Map <pdbId, aaPos> to <transcript, pos>
			//---
			String trid = IdMapperEntry.IDME_TO_REFSEQ.apply(idme);

			// Get transcript
			Transcript tr = trancriptById.get(trid);
			if (tr == null) {
				warn("Transcript not found", trid);
				continue;
			}

			// Transcript's protein doesn't match MSA's protein? Don't use this entry
			if (!checkSequenceGenomeMsas(trid)) continue;

			// Find genomic position based on AA position
			int aa2pos[] = tr.aaNumber2Pos();
			if ((aa2pos.length <= pdb.aaPos) || (pdb.aaPos < 0)) continue;// Position outside amino acid

			// Convert to genomic positions
			int pos = aa2pos[pdb.aaPos];

			//---
			// Map <transcript, pos> to <msaId, aaIdx>
			//---
			// Find MSA and aaIdx
			MsaCoordinates msa = mapTrPos2MsaIdx(trid, pos);
			if (msa != null) {
				// Check sequences
				String seqMsa = msas.getMsa(msa.msaId).getColumnString(msa.msaIdx);

				if ((seqMsa != null) && (pdb.aa == seqMsa.charAt(0))) return msa;
				else {
					Gpr.debug("Error: Mapped sequences do not match: " //
							+ "Seq (MSA): " + seqMsa //
							+ ", Seq (PDB): " + pdb.aa //
							+ "\nPdb: " + pdb //
							+ "\n" + tr //
							+ "\n");
				}
			}
		}
		return null;
	}

	/**
	 * Convert <transcript_id, position> to <msaIdx, aaIdx>
	 */
	public MsaCoordinates mapTrPos2MsaIdx(String trid, int pos) {
		// Find transcript
		Transcript tr = getTranscript(trid);
		if (tr == null) return null;

		// Set genomic coordinates
		GenomicCoordinates gc = new GenomicCoordinates(tr.getChromosome(), pos);

		// Find all MSAs for a given transcript ID
		List<MultipleSequenceAlignment> msaList = getMsas().getMsasByTrId(trid);
		if (msaList == null) return null;

		// Try to map to all MSAs
		for (MultipleSequenceAlignment msa : msaList) {
			// Does this MSA intersect chr:pos?
			if (!msa.intersects(gc)) continue;

			// Try to map to 'msa'
			MsaCoordinates msac = mapMsaTrPos2AaIdx(msa, tr, pos);
			if (msac != null) return msac;
		}

		return null;
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

	public Structure readPdbFile(String pdbFile) {
		try {
			return pdbreader.getStructure(pdbFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
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
