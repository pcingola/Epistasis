package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignmentSet;

/**
 * Find transition matrix by using a maximum likelihood procedure
 *
 * @author pcingola
 */
public class MaxLikelihoodTm {

	public static final int NUM_AA = GprSeq.AMINO_ACIDS.length;
	public static final int NUM_AA_SQUARE = GprSeq.AMINO_ACIDS.length * GprSeq.AMINO_ACIDS.length;

	public static int METHOD = 0;

	boolean verbose = false;
	boolean debug = false;

	int pseudoCount = 1;
	int numSpecies;
	double pi[];
	double time[][];
	LikelihoodTree tree;
	MultipleSequenceAlignmentSet msas;
	TransitionMatrixMarkov Q;
	ArrayRealVector piVect;
	Random random;
	HashMap<String, Double> cacheLogLikelihood;

	public MaxLikelihoodTm(LikelihoodTree tree, MultipleSequenceAlignmentSet msas) {
		this.tree = tree;
		this.msas = msas;
		random = new Random(20140426);
		cacheLogLikelihood = new HashMap<String, Double>();
		numSpecies = tree.childNames().size();
		time = new double[numSpecies][numSpecies];
	}

	/**
	 * Calculate 'stable' probability for each amino acid
	 * Note: We calculate using 'all' alignments and the first alignment
	 */
	public ArrayRealVector calcPi() {
		if (pi != null && piVect != null) return piVect;

		System.out.println("Counting amino acids: ");
		int countAa[] = msas.countAa();
		int tot = 0;
		for (int aa = 0; aa < countAa.length; aa++)
			tot += countAa[aa];

		// Calculate for all AA
		double piAll[];
		piAll = new double[countAa.length];
		if (verbose) System.out.println("AAcode\tAA\tcount_all\tcount_first\tp_all\tp_first");
		for (int aa = 0; aa < countAa.length; aa++) {
			piAll[aa] = countAa[aa] / ((double) tot);
			if (verbose) System.out.println(aa + "\t" + GprSeq.code2aa((byte) aa) + "\t" + countAa[aa] + "\t" + piAll[aa]);
		}

		// Create vector
		pi = piAll;
		piVect = new ArrayRealVector(pi);
		return piVect;
	}

	public RealMatrix calcQ(int seqNum1, int seqNum2) {
		return null;
	}

	/**
	 * Inference of a transition matrix Q
	 */
	public TransitionMatrix estimateTransitionMatrix() {
		calcPi(); // Calculate 'stable' AA distributions

		// Column and row names
		String names[] = new String[GprSeq.AMINO_ACIDS.length];
		for (int i = 0; i < names.length; i++)
			names[i] = GprSeq.code2aa((byte) i) + "";

		//----
		// Estimate matrix
		//---

		// For each pair of species, estimate Q
		System.out.println("Estimate transition matrix");
		int N = GprSeq.AMINO_ACIDS.length;
		Array2DRowRealMatrix QhatSum = new Array2DRowRealMatrix(N, N);
		int count = 0;
		for (int i = 0; i < msas.getNumAligns(); i++) {
			for (int j = i + 1; j < msas.getNumAligns(); j++) {
				TransitionMatrix QhatTmp = estimateTransitionMatrix(i, j, names);
				QhatSum = QhatSum.add(QhatTmp);
				count++;
			}
		}

		// Calculate the average of all estimators
		Q = new TransitionMatrixMarkov(QhatSum.scalarMultiply(1.0 / count));
		Q.setColNames(names);
		Q.setRowNames(names);

		return Q;
	}

	/**
	 * Calculate transition matrix from data
	 * @param seqNum1
	 * @param seqNum2
	 */
	public TransitionMatrix estimateTransitionMatrix(int seqNum1, int seqNum2, String names[]) {
		String seqName1 = msas.getSpecies()[seqNum1];
		String seqName2 = msas.getSpecies()[seqNum2];
		double t = time(seqNum1, seqNum2);
		System.err.println(seqName1 + "\t" + seqName2 + "\ttime: " + t);

		// Count all transitions
		int count[][] = msas.countTransitions(seqNum1, seqNum2);

		// Add pseudo-counts
		for (int i = 0; i < count.length; i++)
			for (int j = 0; j < count.length; j++)
				count[i][j] += pseudoCount;

		// Calculate total counts
		int sum = Arrays.stream(count).flatMapToInt(x -> Arrays.stream(x)).sum();

		// Convert to transition frequencies
		// Estimate matrix P
		double phat[][] = new double[NUM_AA][NUM_AA];
		double n = sum;
		for (int i = 0; i < phat.length; i++)
			for (int j = 0; j < phat.length; j++) {
				// phat[i][j] = (count[i][j] + count[j][i]) / n * pi[i];  // WARNING: Should we be dividing by piAvg here!?!?! Why?

				switch (METHOD) {
				case 0:
					double freq = count[i][j] / n;
					phat[i][j] = freq / pi[i];
					break;

				case 1:
					freq = (count[i][j] + count[j][i]) / (2.0 * n);
					phat[i][j] = freq / pi[i];
					break;

				case 2:
					freq = (count[i][j] + count[j][i]) / (2.0 * n);
					double piAvg = (pi[i] + pi[j]) / 2.0;
					phat[i][j] = freq / piAvg;
					break;
				}
			}

		// Create matrix
		// P(t) = exp(t * Q) = V^T exp(t * D) V  => Q = 1/t log[ P(t) ]
		TransitionMatrixMarkov Phat = new TransitionMatrixMarkov(phat);
		TransitionMatrix Qhat = new TransitionMatrixMarkov(Phat.log().scalarMultiply(1 / t));

		// Remove negative entries from matrix
		double dqhat[][] = Qhat.getData();
		for (int i = 0; i < dqhat.length; i++)
			for (int j = 0; j < dqhat.length; j++)
				if (i != j && dqhat[i][j] < 0) dqhat[i][j] = 0;

		Qhat = new TransitionMatrixMarkov(dqhat);
		Qhat.setColNames(names);
		Qhat.setRowNames(names);
		if (verbose) System.err.println(Gpr.prependEachLine(String.format("\t%s_%s_%.4f\t", seqName1, seqName2, t), Qhat));
		return Qhat;
	}

	public RealVector getPi() {
		return piVect;
	}

	public TransitionMatrix getQ() {
		return Q;
	}

	/**
	 * Load a transition matrix
	 * @param fileName
	 * @return
	 */
	public TransitionMatrix loadTransitionMatrix(String fileName) {
		Q = new TransitionMatrixMarkov(TransitionMatrix.load(fileName));
		return Q;
	}

	/**
	 * Calculate the likelihood (sum all MSAS and pos)
	 */
	public double logLikelyhood() {
		double logLik = 0.0;
		calcPi();

		// Calculate likelihood for each MSA & each base
		for (MultipleSequenceAlignment msa : msas) {
			for (int pos = 0; pos < msa.length(); pos++) {
				if (msa.isSkip(pos)) continue;
				double like = logLikelyhood(msa, pos);
				logLik += -Math.log(like);
			}
			if (verbose) System.out.println("MSA: " + msa.getId() + "\t" + logLik);
		}

		return logLik;
	}

	/**
	 * Calculate log likelihood for an msa:pos
	 * @param msa
	 * @param pos
	 * @return
	 */
	public double logLikelyhood(MultipleSequenceAlignment msa, int pos) {
		// Check cache
		String key = msa.getId() + ":" + pos;
		Double logLik = cacheLogLikelihood.get(key);
		if (logLik != null) return logLik;

		// Set sequence and calculate likelihood
		String seqCol = msa.getColumnString(pos);
		tree.setLeafSequence(seqCol);
		double like = tree.likelihood(Q, pi);
		logLik = -Math.log(like);

		// Cache result
		cacheLogLikelihood.put(key, logLik);

		return logLik;
	}

	/**
	 * Calculate log likelihood for an msa:pos
	 */
	public double logLikelyhood(MultipleSequenceAlignment msa1, int pos1, MultipleSequenceAlignment msa2, int pos2) {
		// Get codes
		byte code1[] = msa1.getColumn(pos1);
		byte code2[] = msa2.getColumn(pos2);

		// Initialize
		int codes[] = new int[code1.length];
		double count[][] = new double[NUM_AA_SQUARE][NUM_AA_SQUARE];
		for (int i = 0; i < code1.length; i++) {
			codes[i] = code1[i] * NUM_AA + code2[i];
			Arrays.fill(count[i], pseudoCount);
		}

		// Count entries
		for (int i = 0; i < code1.length; i++)
			for (int j = i + 1; j < code1.length; j++)
				count[codes[j]][codes[j]] += 1.0;

		// Phat estimation
		// TODO: WRONG!!!
		//		double phat[][] = new double[NUM_AA_SQUARE][NUM_AA_SQUARE];
		//		for (int i = 0; i < count.length; i++)
		//			for (int j = i + 1; j < count.length; j++)
		//				phat[i][j] = (count[i][j] + count[j][i]) / n * pi[i] * pi[i]; // Note: We use symmetry
		//
		//		// Create matrix
		//		// P(t) = exp(t * Q) = V^T exp(t * D) V  => Q = 1/t log[ P(t) ]
		//		TransitionMatrixMarkov Phat = new TransitionMatrixMarkov(phat);
		//		TransitionMatrix Qhat = new TransitionMatrixMarkov(Phat.log().scalarMultiply(1 / t));

		// Calculate likelihood for each lambda

		// Max likelihoood from all lambdas
		return 0;
	}

	/**
	 * Calculate time and cache it
	 */
	double time(int seqNum1, int seqNum2) {
		double t = time[seqNum1][seqNum2];
		if (t > 0 || seqNum1 == seqNum2) return t;
		String seqName1 = msas.getSpecies()[seqNum1];
		String seqName2 = msas.getSpecies()[seqNum2];
		t = time[seqNum1][seqNum2] = tree.distance(seqName1, seqName2);
		return t;
	}
}
