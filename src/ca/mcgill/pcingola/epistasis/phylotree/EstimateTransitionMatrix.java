package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignmentSet;

/**
 * Find transition matrix by using a maximum likelihood procedure
 *
 * @author pcingola
 */
public class EstimateTransitionMatrix {

	public static final int NUM_AA = GprSeq.AMINO_ACIDS.length;
	public static final int NUM_AA_SQUARE = GprSeq.AMINO_ACIDS.length * GprSeq.AMINO_ACIDS.length;

	public static int METHOD = 1;

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
	HashMap<String, Double> cacheLogLikelihood;

	public EstimateTransitionMatrix(LikelihoodTree tree, MultipleSequenceAlignmentSet msas) {
		this.tree = tree;
		this.msas = msas;
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

	/**
	 * Inference of a transition matrix Q
	 */
	public TransitionMatrix estimateTransitionMatrix() {
		calcPi();

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

				// Add all transition matrix estimates
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
					// We normally use this one
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

		// Create matrix
		Qhat = new TransitionMatrixMarkov(dqhat);
		Qhat.setColNames(names);
		Qhat.setRowNames(names);

		// Check
		RealVector z = Qhat.operate(calcPi());
		System.err.println("NORM_QHAT_PI\t" + seqName1 + "\t" + seqName2 + "\t" + t + "\t" + z.getNorm());

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
