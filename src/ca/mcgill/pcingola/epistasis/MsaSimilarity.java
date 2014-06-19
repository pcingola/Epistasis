package ca.mcgill.pcingola.epistasis;

import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Implement a 'similarity' by correlation
 *
 * @author pcingola
 */
public class MsaSimilarity {

	public static final double LOG_2 = Math.log(2.0);

	public static final byte ALIGN_GAP = (byte) -1;
	public static double SHOW_THRESHOLD = 0.99;
	// public static final int MIN_SECOND_TOP_BASE_COUNT = 5;
	public static final int MIN_AA_DISTANCE = 25;
	public static final int SCORE_BINS = 1000;
	public static int SHOW_EVERY = 1000;

	/**
	 * Ratio: number of AA equal the the first / number of non-gap
	 */
	public static double conservation(String seq) {
		if (seq == null || seq.isEmpty()) return 0;

		char fisrt = seq.charAt(0);
		int countEq = 0, count = 0;
		for (int i = 1; i < seq.length(); i++) {
			if (seq.charAt(i) != '-') count++;
			if (seq.charAt(i) == fisrt) countEq++;

		}
		return ((double) countEq) / ((double) count);
	}

	int count = 1;
	protected int minCount = 0;
	protected double threshold = 0;
	protected boolean debug = false;
	protected boolean verbose = true;
	protected int numBases;
	protected double max = 0.0;
	double minScore, maxScore;
	protected int countScore[];
	protected MultipleSequenceAlignmentSet msas;
	protected Random random = new Random();

	public MsaSimilarity(MultipleSequenceAlignmentSet msas) {
		this.msas = msas;
		minScore = 0.0;
		maxScore = 1.0;
		countScore = new int[SCORE_BINS];
		numBases = 1;
	}

	/**
	 * Calculate one sample of a random distribution
	 */
	void backgroundDistribution() {
		while (true) {
			// Select one MSA and position randomly
			MultipleSequenceAlignment msai = msas.rand(random);
			int posi = msai.randomColumnNumber(random);
			if (msai.isSkip(posi)) continue;

			// Select another MSA and position randomly
			String trId = msai.getTranscriptId();
			List<MultipleSequenceAlignment> msasTr = msas.getMsas(trId);
			MultipleSequenceAlignment msaj = msasTr.get(random.nextInt(msasTr.size()));
			int posj = msaj.randomColumnNumber(random);
			if (msaj.isSkip(posj)) continue;

			// Same MSA and position? Find another random
			if (posi == posj && msai.getId().equals(msaj.getId())) continue;

			// Calculate
			double calc = calc(msai, msaj, posi, posj);
			if (debug) System.err.println(calc + "\t" + showSeqs(msai, msaj, posi, posj));
			else if (verbose) System.out.println(calc);
			else Gpr.showMark(count++, SHOW_EVERY);
			return;
		}
	}

	/**
	 * Measure similarity between all alignments
	 */
	public void backgroundDistribution(int numberOfSamples) {
		// Pre-calculate skip on all msas
		msas.calcSkip();

		// Calculate in parallel
		Timer.showStdErr("Calculating " + numberOfSamples + " iterations");
		IntStream.range(1, numberOfSamples) //
				.parallel() //
				.forEach(i -> backgroundDistribution());
	}

	/**
	 * Measure similarity: Correlation between two loci
	 */
	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		int count = 0, sum = 0;
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
				matchBases[i] = basei;
			}
		}

		// Only a few organisms align? Filter out
		if (count < minCount) return Double.NaN;

		// Results
		double corr = ((double) sum) / ((double) count);

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

	public void setThreshold(double threshold) {
		this.threshold = threshold;
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
		sb.append(msai.getId() + " [ " + posi + " ]\t" + msaj.getId() + " [ " + posj + " ]");

		int numAligns = msai.getNumSeqs();
		for (int n = 0; n < numBases; n++) {
			sb.append("\t" + n + "\t");
			for (int i = 0; i < numAligns; i++)
				sb.append(msai.getChar(i, posi + n));
		}

		numAligns = msaj.getNumSeqs();
		for (int n = 0; n < numBases; n++) {
			sb.append("\t" + n + "\t");
			for (int i = 0; i < numAligns; i++)
				sb.append(msaj.getChar(i, posj + n));
		}

		return sb.toString();
	}

	//	void similarity(MultipleSequenceAlignment msai, List<MultipleSequenceAlignment> msasTr) {
	//		msasTr.stream() //
	//				.filter(msaj -> msai.compareTo(msaj) <= 0) //
	//				.forEach(msaj -> similarity(msai, msaj)) //
	//		;
	//	}
	//
	//	/**
	//	 * Measure correlation between two alignments
	//	 */
	//	void similarity(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj) {
	//		int maxi = msai.getSeqLen();
	//		int maxj = msaj.getSeqLen();
	//
	//		for (int posi = 0; posi < maxi; posi++) {
	//			if (!msai.isSkip(posi)) {
	//				int posjmin = 0;
	//
	//				// Same alignment? Then make sure we are at least MIN_AA_DISTANCE away (obviously two adjacent AAs are correlated and interacting)
	//				if (msai.getId().equals(msaj.getId())) posjmin = posi + MIN_AA_DISTANCE;
	//
	//				for (int posj = posjmin; posj < maxj; posj++)
	//					if (!msaj.isSkip(posj)) {
	//						try {
	//							similarity(msai, msaj, posi, posj);
	//						} catch (Throwable t) {
	//							// Show error details
	//							Gpr.debug("ERROR processing:" //
	//									+ "\n\tmsa i : " + msai.getId() //
	//									+ "\n\tpos i : " + posi //
	//									+ "\n\tlen i : " + msai.getSeqLen() //
	//									+ "\n" //
	//									+ "\n\tmsa j : " + msaj.getId() //
	//									+ "\n\tpos j : " + posj //
	//									+ "\n\tlen j : " + msaj.getSeqLen() //
	//							);
	//							throw new RuntimeException(t);
	//						}
	//					}
	//			}
	//		}
	//	}

	//	/**
	//	 * Measure correlation between two loci
	//	 */
	//	protected void similarity(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
	//		double calc = calc(msai, msaj, posi, posj);
	//		if (Double.isNaN(calc)) return;
	//		if (debug) System.err.println(calc + "\t" + showSeqs(msai, msaj, posi, posj));
	//		else Gpr.showMark(count++, SHOW_EVERY);
	//	}

	//	/**
	//	 * Correlation with another multiple sequence alignment
	//	 * @param msai
	//	 */
	//	public void similarity(String trId) {
	//		// Get all MSAs for this transcript ID
	//		List<MultipleSequenceAlignment> msasTr = msas.getMsas(trId);
	//		msasTr.stream().forEach(msa -> similarity(msa, msasTr));
	//	}

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
