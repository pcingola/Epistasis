package ca.mcgill.pcingola.epistasis;

import java.util.Random;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Represents an alignment of multiple sequnces of Amini Acids
 *
 * Note: This class does NOT perform an alignmet, only represents the result loaded form a file
 *
 * @author pcingola
 */
public class MultipleSequenceAlignment implements Comparable<MultipleSequenceAlignment> {

	public static final double MAX_GAP_PERCENT = 0.3;
	public static final int ROTATE_BITS = 5;

	String transcriptId;
	String chromo;
	int start, end;
	boolean strandNegative;
	byte[][] align;
	Boolean skip[];

	public MultipleSequenceAlignment(String id, int numAlign, int length) {
		transcriptId = id;
		align = new byte[numAlign][length];
	}

	public synchronized void calcSkip() {
		for (int i = 0; i < getSeqLen(); i++)
			isSkip(i);
	}

	@Override
	public int compareTo(MultipleSequenceAlignment msa) {
		int comp = Chromosome.compare(chromo, msa.chromo);
		if (comp > 0) return Integer.MAX_VALUE;
		if (comp < 0) return Integer.MIN_VALUE;

		if (start > msa.end) return start - msa.end;
		if (end < msa.start) return end - msa.start;

		return 0;
	}

	/**
	 * Count amino acids for a specific alignment
	 * @param alignNum
	 * @param countAa
	 */
	public void countAa(int alignNum, int[] countAa) {
		int len = length();
		for (int j = 0; j < len; j++) {
			byte base = align[alignNum][j];
			if (base >= 0) countAa[base]++;
		}
	}

	/**
	 * Count amino acids in all sequences
	 * @param countAa : vector to be updated
	 * @return
	 */
	public void countAa(int[] countAa) {
		for (int i = 0; i < align.length; i++)
			countAa(i, countAa);
	}

	/**
	 * How many 'changes' happen in this position?
	 * @param pos
	 * @return
	 */
	public int countChanges(int pos) {
		int changes = 0;
		byte prevBase = ' ';

		for (int i = 0; i < getNumSeqs(); i++) {
			byte base = getCode(i, pos);
			if (base < 0) continue;

			if (prevBase != base) changes++;
			prevBase = base;
		}

		return changes - 1;
	}

	/**
	 * Count number of transitions between two sequences
	 * @param seqNum1
	 * @param seqNum2
	 * @return
	 */
	public void countTransitions(int seqNum1, int seqNum2, int counts[][]) {
		int len = length();
		for (int j = 0; j < len; j++) {
			byte base1 = align[seqNum1][j];
			if (base1 < 0) continue;

			byte base2 = align[seqNum2][j];
			if (base2 < 0) continue;

			counts[base1][base2]++;
		}
	}

	/**
	 * Percentage of gaps at position 'pos'
	 * @param pos
	 * @return
	 */
	double gapPercent(int pos) {
		int gaps = 0;
		for (int i = 0; i < getNumSeqs(); i++) {
			byte base = getCode(i, pos);
			if (base < 0) gaps++;
		}
		return ((double) gaps) / ((double) getNumSeqs());
	}

	public byte[][] getAlign() {
		return align;
	}

	public char getChar(int sequenceNum, int baseNum) {
		return GprSeq.code2aa(align[sequenceNum][baseNum]);
	}

	public String getChromo() {
		return chromo;
	}

	public byte getCode(int sequenceNum, int baseNum) {
		return align[sequenceNum][baseNum];
	}

	public int getCodeInt(int sequenceNum, int baseStart, int baseCount) {
		int tot = 0;
		for (int i = 0; i < baseCount; i++) {
			byte base = align[sequenceNum][baseStart + i];
			if (base == MsaSimilarity.ALIGN_GAP) return 0;
			tot = (tot << ROTATE_BITS) | base;
		}

		return tot;
	}

	public long getCodeLong(int sequenceNum, int baseStart, int baseCount) {
		long tot = 0;
		for (int i = 0; i < baseCount; i++)
			tot = (tot << 5) | align[sequenceNum][baseStart + i];
		return tot;
	}

	/**
	 * Get all characters in column 'colNum'
	 * @param colNum
	 * @return
	 */
	public byte[] getColumn(int colNum) {
		byte col[] = new byte[size()];
		for (int i = 0; i < col.length; i++)
			col[i] = align[i][colNum];

		return col;
	}

	/**
	 * Get all characters in column 'colNum'
	 */
	public String getColumnString(int colNum) {
		if (colNum < 0 || colNum >= length()) return null;

		char col[] = new char[size()];
		for (int i = 0; i < col.length; i++)
			col[i] = getChar(i, colNum);

		return new String(col);
	}

	public int getEnd() {
		return end;
	}

	public String getId() {
		return transcriptId + "_" + chromo + ":" + start + "-" + end;
	}

	/**
	 * Number of sequences (i.e. number of species aligned)
	 */
	public int getNumSeqs() {
		return align.length;
	}

	/**
	 * Get all characters in row 'rowNum'
	 */
	public String getRowString(int rowNum) {
		if (rowNum < 0 || rowNum >= getNumSeqs()) return null;

		char row[] = new char[length()];
		for (int i = 0; i < row.length; i++)
			row[i] = getChar(rowNum, i);

		return new String(row);
	}

	/**
	 * Sequence length (i.e. number of columns)
	 */
	public int getSeqLen() {
		return align[0].length;
	}

	public int getStart() {
		return start;
	}

	public String getTranscriptId() {
		return transcriptId;
	}

	/**
	 * Is position 'pos' fully conserved?
	 */
	public boolean isFullyConserved(int pos) {
		byte prevBase = -1;
		for (int i = 0; i < getNumSeqs(); i++) {
			byte base = getCode(i, pos);
			if (base < 0) continue;
			if (prevBase != base && prevBase >= 0) return false;
			prevBase = base;
		}

		return true;
	}

	/**
	 * Are the 'n' bases at position 'pos' fully conserved?
	 */
	public boolean isFullyConserved(int pos, int n) {
		for (int i = pos - n; i < (pos + n); i++)
			if ((i >= 0) && !isFullyConserved(i)) return false;

		return true;
	}

	/**
	 * Should we skip position 'pos'?
	 * @param pos
	 * @return
	 */
	public synchronized boolean isSkip(int pos) {
		if (skip == null) skip = new Boolean[getSeqLen()];
		if (skip[pos] != null) return skip[pos];

		// Should we skip this base?
		boolean skipThis = false;

		if (getCode(0, pos) < 0) skipThis = true; // Is it a GAP in the first alignment (human). Skip, because we only care about human
		if (!skipThis) skipThis = (gapPercent(pos) >= MAX_GAP_PERCENT);
		if (!skipThis) skipThis = isFullyConserved(pos);
		skip[pos] = skipThis;

		return skip[pos];
	}

	public boolean isStrandNegative() {
		return strandNegative;
	}

	public boolean isStrandPositive() {
		return !strandNegative;
	}

	/**
	 * Length = Sequence length (i.e. number of columns)
	 */
	public int length() {
		return getSeqLen();
	}

	/**
	 * Return a random column number
	 */
	public int randomColumnNumber(Random random) {
		return random.nextInt(getSeqLen());
	}

	/**
	 * Set sequence number 'num'
	 * @param seqNum
	 * @param seq
	 */
	public void set(int seqNum, String seq) {
		for (int i = 0; i < seq.length(); i++)
			align[seqNum][i] = GprSeq.aa2Code(seq.charAt(i));
	}

	/**
	 * Set genomic coordinates (chr:start-end)
	 * @param chr
	 * @param start
	 * @param end
	 */
	public void set(String chr, int start, int end, boolean strandNegative) {
		chromo = Chromosome.simpleName(chr);
		this.start = start;
		this.end = end;
		this.strandNegative = strandNegative;
	}

	/**
	 * Size = Number of sequences (i.e. number of species aligned)
	 */
	public int size() {
		return getNumSeqs();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getId() + "\n");

		for (int i = 0; i < align.length; i++) {
			sb.append("\t");
			for (int j = 0; j < align[i].length; j++)
				sb.append(GprSeq.code2aa(align[i][j]));
			sb.append("\t" + i + "\n");
		}

		return sb.toString();
	}
}
