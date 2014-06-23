package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Load a multiple sequence alignment file (UCSC)
 * E.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/alignments/knownCanonical.exonAA.fa.gz
 *
 * @author pcingola
 */
public class MultipleSequenceAlignmentSet implements Iterable<MultipleSequenceAlignment> {

	public static boolean debug = false;
	public static boolean verbose = false;
	public double SHOW_THRESHOLD = 0.99;
	public final int MIN_COUNT_THRESHOLD = 50;
	public final int MIN_SECOND_TOP_BASE_COUNT = 5;

	ArrayList<MultipleSequenceAlignment> msas;
	AutoHashMap<String, List<MultipleSequenceAlignment>> msasByTrId;
	HashMap<String, MultipleSequenceAlignment> msaById;
	int numAligns;
	String sequenceAlignmentFile;
	String species[];

	public MultipleSequenceAlignmentSet(String sequenceAlignmentFile, int numAligns) {
		this.numAligns = numAligns;
		this.sequenceAlignmentFile = sequenceAlignmentFile;
		species = new String[numAligns];
		msas = new ArrayList<MultipleSequenceAlignment>();
		msasByTrId = new AutoHashMap<String, List<MultipleSequenceAlignment>>(new ArrayList<MultipleSequenceAlignment>());
		msaById = new HashMap<String, MultipleSequenceAlignment>();
	}

	public void calcSkip() {
		Timer.showStdErr("Pre-calculating skips.");
		getMsas().parallelStream().forEach(MultipleSequenceAlignment::calcSkip);
	}

	/**
	 * Count number of amino acids
	 * @param countAa
	 */
	public int[] countAa() {
		int counts[] = new int[GprSeq.AMINO_ACIDS.length];
		forEach(m -> m.countAa(counts));
		return counts;
	}

	/**
	 * Count number of amino acids for a specific alignment
	 * @param countAa
	 */
	public int[] countAa(int alignNum) {
		int counts[] = new int[GprSeq.AMINO_ACIDS.length];
		forEach(m -> m.countAa(alignNum, counts));
		return counts;
	}

	/**
	 * Count number of transitions between two sequences
	 */
	public int[][] countTransitions(int seqNum1, int seqNum2) {
		int counts[][] = new int[GprSeq.AMINO_ACIDS.length][GprSeq.AMINO_ACIDS.length];
		forEach(m -> m.countTransitions(seqNum1, seqNum2, counts));
		return counts;
	}

	/**
	 * Return '2 * numBases + 1' string representing the column sequences
	 * at position msaId:pos and the surrounding 'numBases' columns
	 */
	public String[] findColSequences(String msaId, int pos, int numBases) {
		// Initialize
		String seqs[] = new String[2 * numBases + 1];

		// Find positions 'pos' and after
		MultipleSequenceAlignment msa = getMsa(msaId);
		int maxj = 2 * numBases + 1;
		for (int i = pos, j = numBases; j < maxj; i++, j++) {
			if (i >= msa.length()) {
				msa = msaTranscriptNext(msa);
				if (msa == null) break;
				i = 0;
			}
			seqs[j] = msa.getColumnString(i);
		}

		// Find positions before 'pos'
		msa = getMsa(msaId);
		for (int i = pos - 1, j = numBases - 1; j >= 0; i--, j--) {
			if (i < 0) {
				msa = msaTranscriptPrev(msa);
				if (msa == null) break;
				i = msa.length() - 1;
			}
			seqs[j] = msa.getColumnString(i);
		}

		return seqs;
	}

	/**
	 * Find a row sequence
	 * @param trid : Transcript ID
	 * @param chr : Chromosome name (to check that matches the MSA). It can be null, in which case, checking is skipped.
	 * @return
	 */
	public String findRowSequence(String trid, String chr) {
		// Find all MSA
		List<MultipleSequenceAlignment> msaList = msasByTrId.get(trid);
		if (msaList == null) return null;

		// Get all msas for this 'trid'
		ArrayList<MultipleSequenceAlignment> msasTr = new ArrayList<>();
		boolean reverse = false;
		for (MultipleSequenceAlignment msa : msaList) {
			// Different chromosome or position? Skip
			if (chr != null && !msa.getChromo().equals(chr)) continue;
			msasTr.add(msa);
			reverse |= msa.isStrandNegative();
		}

		// Sort
		if (reverse) Collections.sort(msasTr, Collections.reverseOrder());
		else Collections.sort(msasTr);

		StringBuilder sb = new StringBuilder();
		for (MultipleSequenceAlignment msa : msasTr)
			sb.append(msa.getRowString(0));

		return sb.toString();
	}

	public MultipleSequenceAlignment getMsa(String msaId) {
		return msaById.get(msaId);
	}

	public ArrayList<MultipleSequenceAlignment> getMsas() {
		return msas;
	}

	public List<MultipleSequenceAlignment> getMsas(String trId) {
		return msasByTrId.get(trId);
	}

	public int getNumAligns() {
		return numAligns;
	}

	public String[] getSpecies() {
		return species;
	}

	@Override
	public Iterator<MultipleSequenceAlignment> iterator() {
		return msas.iterator();
	}

	/**
	 * @param args
	 */
	public void load() {
		Timer.showStdErr("Loading file '" + sequenceAlignmentFile + "'");
		LineFileIterator lif = new LineFileIterator(sequenceAlignmentFile);

		while (lif.hasNext()) {
			int seqLen = -1;
			MultipleSequenceAlignment msa = null;

			// Read an alignment of a protein
			for (int i = 0; i < numAligns; i++) {
				String header = lif.next();

				// There might be an extra empty line in the first header
				while (header != null && header.isEmpty())
					header = lif.next();
				if (header == null) break;

				// Is this a header?
				if (!header.startsWith(">")) throw new RuntimeException("Error (file '" + sequenceAlignmentFile + "', line " + lif.getLineNum() + "): Expecting header empty line, got: '" + header + "'");

				// Parse and check species
				String fields[] = header.split("_");
				String transcriptId = fields[0].substring(1) + "_" + fields[1];
				String speciesName = fields[2];
				if (species[i] == null) species[i] = speciesName;
				else if (!speciesName.equals(species[i])) throw new RuntimeException("Error (file '" + sequenceAlignmentFile + "', line " + lif.getLineNum() + "): Expecting species '" + species[i] + "', got: '" + speciesName + "'");

				// Read sequence
				String sequence = lif.next();
				if (sequence == null) break;

				// Check sequence length
				if (seqLen < 0) {
					seqLen = sequence.length();

					// Chr:start=end
					fields = header.split(" ");
					String chrpos = fields[4];
					int idxPos = chrpos.indexOf(':');
					int idxEnd = chrpos.indexOf('-');
					String chr = chrpos.substring(0, idxPos);
					String posStart = chrpos.substring(idxPos + 1, idxEnd);
					String posEnd = chrpos.substring(idxEnd + 1, chrpos.length() - 1);
					boolean strand = (chrpos.charAt(chrpos.length() - 1) == '-');
					int start = Gpr.parseIntSafe(posStart) - 1;
					int end = Gpr.parseIntSafe(posEnd) - 1;

					if (debug) System.out.println(transcriptId + " " + chr + ":" + start + "-" + end);
					msa = new MultipleSequenceAlignment(transcriptId, numAligns, seqLen);
					msa.set(chr, start, end, strand);
				} else if (sequence.length() != seqLen) throw new RuntimeException("Error (file '" + sequenceAlignmentFile + "', line " + lif.getLineNum() + "): Expecting sequence of length " + seqLen);

				// Add sequence
				if (debug) System.out.println("\t" + sequence + "\t" + species[i]);

				// Remove ambiguous amino acids
				sequence = sequence.replace('B', '-');
				sequence = sequence.replace('Z', '-');
				sequence = sequence.replace('J', '-');
				sequence = sequence.replace('X', '-');

				// Remove rare amino acids
				sequence = sequence.replace('U', '-');
				sequence = sequence.replace('O', '-');

				// Stop codons
				sequence = sequence.replace('*', '-');

				// Set sequence
				msa.set(i, sequence);
			}

			if (msa != null) {
				msas.add(msa);
				msasByTrId.getOrCreate(msa.getTranscriptId()).add(msa);
				msaById.put(msa.getId(), msa);
				if (verbose) System.out.println(msa.getId());
			}

			// Empty line separator
			String emptyLine = lif.next();
			if (emptyLine != null && !emptyLine.isEmpty()) throw new RuntimeException("Error (file '" + sequenceAlignmentFile + "', line " + lif.getLineNum() + "): Expecting an empty line!");
		}

		// Sort lists
		sortTranscriptLists();

		Timer.showStdErr("Done. Total number of alignments: " + msas.size());
	}

	/**
	 * Obtain next MSA (for the same transcript)
	 */
	MultipleSequenceAlignment msaTranscriptNext(MultipleSequenceAlignment msa) {
		return null;
	}

	/**
	 * Obtain previous MSA (for the same transcript)
	 */
	MultipleSequenceAlignment msaTranscriptPrev(MultipleSequenceAlignment msa) {
		return null;
	}

	/**
	 * Get a random element
	 */
	public MultipleSequenceAlignment rand(Random random) {
		return msas.get(random.nextInt(msas.size()));
	}

	public int size() {
		return msas.size();
	}

	/**
	 * Get all transcript lists sorted by strand
	 */
	void sortTranscriptLists() {
		for (List<MultipleSequenceAlignment> l : msasByTrId.values()) {
			if (l.isEmpty()) continue;
			MultipleSequenceAlignment msa = l.get(0); // First msa in the list

			// Sort by strand
			if (msa.isStrandPositive()) Collections.sort(l);
			else Collections.sort(l, Comparator.reverseOrder());
		}
	}

}
