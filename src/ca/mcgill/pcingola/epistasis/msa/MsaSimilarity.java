package ca.mcgill.pcingola.epistasis.msa;

import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.entropy.EntropySeq;
import ca.mcgill.pcingola.epistasis.entropy.EntropySeq.InformationFunction;
import ca.mcgill.pcingola.epistasis.pdb.DistanceResult;

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
	protected InformationFunction func = InformationFunction.MI;

	public MsaSimilarity(MultipleSequenceAlignmentSet msas, InformationFunction func) {
		this.msas = msas;
		this.func = func;
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
			List<MultipleSequenceAlignment> msasTr = msas.getMsasByTrId(trId);
			MultipleSequenceAlignment msaj = msasTr.get(random.nextInt(msasTr.size()));
			int posj = msaj.randomColumnNumber(random);
			if (msaj.isSkip(posj)) continue;

			// Same MSA and position? Find another random
			if (posi == posj && msai.getId().equals(msaj.getId())) continue;

			// Calculate
			double calc = calc(msai, msaj, posi, posj);
			if (debug) System.err.println(calc + "\t" + showSeqs(msai, msaj, posi, posj));
			else if (verbose) System.out.printf("%.6e\n", calc);
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

	public double calc(DistanceResult d) {
		MultipleSequenceAlignment msai = msas.getMsa(d.msa1);
		if (msai == null) return Double.NaN;

		MultipleSequenceAlignment msaj = msas.getMsa(d.msa2);
		if (msaj == null) return Double.NaN;

		return calc(msai, msaj, d.msaIdx1, d.msaIdx2);
	}

	/**
	 * Measure similarity: Correlation between two loci
	 */
	public double calc(MultipleSequenceAlignment msai, MultipleSequenceAlignment msaj, int posi, int posj) {
		int numAligns = msas.getNumAligns();
		byte basesi[] = new byte[numAligns];
		byte basesj[] = new byte[numAligns];

		// Convert
		for (int i = 0; i < numAligns; i++) {
			basesi[i] = msai.getCode(i, posi);
			basesj[i] = msaj.getCode(i, posj);
		}

		// Calculate
		double score = 0;
		switch (func) {
		case MI:
			score = EntropySeq.mutualInformation(basesi, basesj);
			break;

		case VARINF:
			score = EntropySeq.variationOfInformation(basesi, basesj);
			break;

		case HXY:
			score = EntropySeq.entropy(basesi, basesj);
			break;

		default:
			throw new RuntimeException("Unknown function type " + func);
		}

		incScore(score);
		return score;
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
