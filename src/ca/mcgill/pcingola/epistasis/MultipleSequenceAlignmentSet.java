package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
	public static byte ALIGN_GAP = (byte) 0;
	public double SHOW_THRESHOLD = 0.99;
	public final int MIN_COUNT_THRESHOLD = 50;
	public final int MIN_SECOND_TOP_BASE_COUNT = 5;

	ArrayList<MultipleSequenceAlignment> msas;
	AutoHashMap<String, List<MultipleSequenceAlignment>> msasById;
	int numAligns;
	String sequenceAlignmentFile;
	String species[];

	public MultipleSequenceAlignmentSet(String sequenceAlignmentFile, int numAligns) {
		this.numAligns = numAligns;
		this.sequenceAlignmentFile = sequenceAlignmentFile;
		species = new String[numAligns];
		msas = new ArrayList<MultipleSequenceAlignment>();
		msasById = new AutoHashMap<String, List<MultipleSequenceAlignment>>(new ArrayList<MultipleSequenceAlignment>());
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
	 * @param seqNum1
	 * @param seqNum2
	 * @return
	 */
	public int[][] countTransitions(int seqNum1, int seqNum2) {
		int counts[][] = new int[GprSeq.AMINO_ACIDS.length][GprSeq.AMINO_ACIDS.length];
		forEach(m -> m.countTransitions(seqNum1, seqNum2, counts));
		return counts;
	}

	/**
	 * Find a multiple sequence alignment
	 */
	public String findColumnSequence(String trid, String chr, int pos) {
		// Find all MSA
		List<MultipleSequenceAlignment> msaList = msasById.get(trid);
		if (msaList == null) return null;

		// Check all MSA
		for (MultipleSequenceAlignment msa : msaList) {
			// Different chromosome?
			if (!msa.getChromo().equals(chr)) {
				Gpr.debug("Transcript '" + trid + "' is expected in chromosome '" + chr + "', but found in chromosome '" + msa.getChromo() + "'. Skipped");
				continue;
			}

			// Different position?
			if (pos < msa.getStart() || msa.getEnd() < pos) {
				Gpr.debug("Transcript '" + trid + "' is expected in position " + pos + ", but found in interval [ " + msa.getStart() + " , " + msa.getEnd() + " ]. Skipped");
				continue;
			}

			// Found!
			int idx = pos - msa.getStart();
			return msa.getColumnString(idx);
		}

		return null;
	}

	public ArrayList<MultipleSequenceAlignment> getMsas() {
		return msas;
	}

	public List<MultipleSequenceAlignment> getMsas(String trid) {
		return msasById.get(trid);
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
					int start = Gpr.parseIntSafe(posStart);
					int end = Gpr.parseIntSafe(posEnd);

					if (debug) System.out.println(transcriptId + " " + chr + ":" + start + "-" + end);
					msa = new MultipleSequenceAlignment(transcriptId, numAligns, seqLen);
					msa.set(chr, start, end);
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
				msasById.getOrCreate(msa.getTranscriptId()).add(msa);
				if (verbose) System.out.println(msa.getId());
			}

			// Empty line separator
			String emptyLine = lif.next();
			if (emptyLine != null && !emptyLine.isEmpty()) throw new RuntimeException("Error (file '" + sequenceAlignmentFile + "', line " + lif.getLineNum() + "): Expecting an empty line!");
		}

		Timer.showStdErr("Done. Total number of alignments: " + msas.size());
	}

	public int size() {
		return msas.size();
	}

}
