package ca.mcgill.pcingola.epistasis;

import java.util.stream.Collectors;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Implement a 'similarity' by correlation
 *
 * @author pcingola
 */
public class MsaSimilarity {

	public static final byte ALIGN_GAP = (byte) -1;
	public static double SHOW_THRESHOLD = 0.99;
	public static final int MIN_COUNT_THRESHOLD = 50;
	public static final int MIN_SECOND_TOP_BASE_COUNT = 5;
	public static final int MIN_AA_DISTANCE = 10;
	public static final int SCORE_BINS = 1000;

	protected boolean debug = false;
	protected int numBases;
	protected double max = 0.0;
	double minScore, maxScore;
	protected int countScore[];
	protected MultipleSequenceAlignmentSet msas;

	public MsaSimilarity(MultipleSequenceAlignmentSet msas) {
		this.msas = msas;
		minScore = 0.0;
		maxScore = 1.0;
		countScore = new int[SCORE_BINS];
		numBases = 1;
	}

	/**
	 * Measure similarity: Correlation between two loci
	 */
	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		int count = 0, sum = 0;
		int changes = 0;
		byte prevBase = -1;
		int numAligns = msas.getNumAligns();
		byte matchBases[] = new byte[numAligns];

		// Count matching bases
		for (int i = 0; i < numAligns; i++) {
			byte basei = msai.getCode(i, posi);
			byte basej = msaj.getCode(i, posj);

			// TODO: Take gaps into account?
			if ((basei == ALIGN_GAP) || (basej == ALIGN_GAP)) {
				matchBases[i] = '-';
				continue;
			}

			count++;
			if (basei == basej) {
				sum++;
				if (prevBase != basei) changes++;
				matchBases[i] = prevBase = basei;
			}
		}

		// Only a few organisms align? Filter out
		if (count < MIN_COUNT_THRESHOLD) return Double.NaN;

		// Filter out if they match only in the same base
		if (changes <= 1) return Double.NaN;

		// Filter out is too low count on the second 'top' base
		if (secondMostCommonBaseCount(matchBases) < MIN_SECOND_TOP_BASE_COUNT) return Double.NaN;

		// Results
		float corr = ((float) sum) / ((float) count);

		incScore(corr);
		if (corr > SHOW_THRESHOLD) return corr;
		return Double.NaN;
	}

	/**
	 * Increment counter for score 'score'
	 * @param score
	 */
	protected void incScore(double score) {
		int idx = mapScore(score);
		if (idx >= countScore.length) throw new RuntimeException("score:" + score + " out of range. maxScore: " + maxScore);
		countScore[idx]++;
	}

	protected boolean isRecord(double score) {
		if (score > max) {
			max = score;
			return true;
		}
		return false;
	}

	protected int mapScore(double score) {
		return (int) (SCORE_BINS * ((score - minScore) / (maxScore - minScore)));
	}

	/**
	 * What is the count for the second most common base
	 * @param pos
	 * @return
	 */
	public int secondMostCommonBaseCount(byte bases[]) {
		int count[] = new int[256];
		for (int i = 0; i < count.length; i++)
			count[i] = 0;

		for (int i = 0; i < bases.length; i++)
			count[bases[i]]++;

		int max = 0, second = 0;
		for (int i = 0; i < count.length; i++) {
			if (i == ' ' || i == MsaSimilarity.ALIGN_GAP) continue;

			if (count[i] > max) max = count[i];
			else if (count[i] > second) second = count[i];
		}

		return second;
	}

	/**
	 * Show columns from MultipleSequenceAlignment
	 * @param msai
	 * @param msaj
	 * @param posi
	 * @param posj
	 * @return
	 */
	public String showSeqs(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		StringBuilder sb = new StringBuilder();
		sb.append(msai.getId() + " [ " + posi + " ]\t" + msaj.getId() + " [ " + posj + " ]\n");

		int numAligns = msai.getNumSeqs();
		for (int n = 0; n < numBases; n++) {
			sb.append("\n\t" + n + "\t");
			for (int i = 0; i < numAligns; i++)
				sb.append(msai.getChar(i, posi + n));
		}

		sb.append("\n");

		numAligns = msaj.getNumSeqs();
		for (int n = 0; n < numBases; n++) {
			sb.append("\n\t" + n + "\t");
			for (int i = 0; i < numAligns; i++)
				sb.append(msaj.getChar(i, posj + n));
		}
		sb.append('\n');

		return sb.toString();
	}

	/**
	 * Measure similarity between all alignments
	 */
	public void similarity() {
		// Pre-calculate skip on all msas
		Timer.showStdErr("Pre-calculating skips.");
		msas.getMsas().parallelStream().forEach(MultipleSequenceAlignment::calcSkip);
		Timer.showStdErr("Done.");

		// Correlation
		Timer.showStdErr("Calculating correlations.");
		msas.getMsas().parallelStream().map(msa -> similarity(msa)).forEach(System.out::print);
		Timer.showStdErr("Done.");
	}

	/**
	 * Correlation with another multiple sequence alignment
	 * @param msai
	 */
	public String similarity(MultipleSequenceAlignment msai) {
		String result = msas.getMsas().stream() //
				.filter(msaj -> msai.getTranscriptId().equals(msaj.getTranscriptId())) // Filter by name (only same sequence)
				.filter(msaj -> msai.getId().compareTo(msaj.getId()) <= 0) // Filter by name (don't work twice on the same pair of sequences)
				.map(msaj -> similarity(msai, msaj)) // Calculate similarities
				.collect(Collectors.joining()) // Join results to one string
				;

		return result.isEmpty() ? "" : msai.getId() + "\n" + result;
	}

	/**
	 * Measure correlation between two alignments
	 */
	public String similarity(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj) {
		int maxi = msai.getSeqLen();
		int maxj = msaj.getSeqLen();

		StringBuilder sb = new StringBuilder();
		for (int posi = 0; posi < maxi; posi++) {
			if (!msai.isSkip(posi)) {
				int posjmin = 0;

				// Same alignment? Then make sure we are at least MIN_AA_DISTANCE away (obviously two adjacent AAs are correlated and interacting)
				if (msai.getId().equals(msaj.getId())) posjmin = posi + MIN_AA_DISTANCE;

				for (int posj = posjmin; posj < maxj; posj++)
					if (!msaj.isSkip(posj)) {
						try {
							sb.append(similarity(msai, msaj, posi, posj));
						} catch (Throwable t) {
							// Show error details
							Gpr.debug("ERROR processing:" //
									+ "\n\tmsa i : " + msai.getId() //
									+ "\n\tpos i : " + posi //
									+ "\n\tlen i : " + msai.getSeqLen() //
									+ "\n" //
									+ "\n\tmsa j : " + msaj.getId() //
									+ "\n\tpos j : " + posj //
									+ "\n\tlen j : " + msaj.getSeqLen() //
									);
							throw new RuntimeException(t);
						}
					}
			}
		}

		return sb.toString();
	}

	/**
	 * Measure correlation between two loci
	 */
	public String similarity(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		double calc = calc(msai, msaj, posi, posj);
		if (Double.isNaN(calc)) return "";

		String message = "\t" + calc + "\t" + showSeqs(msai, msaj, posi, posj);

		if (debug) System.err.println(message);

		// Show if record
		if (isRecord(calc)) System.err.println("New record: " + calc + "\n" + message);

		return message.toString();
	}

	public int size() {
		return msas.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("score\tcount\n");
		for (int i = 0; i < countScore.length; i++) {
			double sc = minScore + ((maxScore - minScore) / SCORE_BINS) * i;
			sb.append(sc + "\t" + countScore[i] + "\n");
		}

		return sb.toString();
	}
}
