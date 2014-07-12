package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ca.mcgill.mcb.pcingola.stats.Counter;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Tuple;
import ca.mcgill.pcingola.epistasis.MultipleSequenceAlignmentSet;

/**
 * Find transition matrix by using a maximum likelihood procedure
 *
 * @author pcingola
 */
public class EstimateTransitionMatrix {

	public static final int NUM_AA = GprSeq.AMINO_ACIDS.length;
	public static final int NUM_AA_SQUARE = GprSeq.AMINO_ACIDS.length * GprSeq.AMINO_ACIDS.length;

	public static int METHOD = 0;
	public static int REMOVE_NEGATIVES = 0;
	public static int PSEUDO_COUNTS = 0;

	boolean verbose = false;
	boolean debug = false;
	int N;
	int pseudoCount = 1;
	int numSpecies;
	double pi[];
	double time[][];
	String names[];
	LikelihoodTreeAa tree;
	MultipleSequenceAlignmentSet msas;
	TransitionMatrixMarkov Q;
	ArrayRealVector piVect;
	HashMap<String, Double> cacheLogLikelihood;

	public static String methods() {
		return "METHOD:" + EstimateTransitionMatrix.METHOD //
				+ "_RMNEGS:" + EstimateTransitionMatrix.REMOVE_NEGATIVES //
				+ "_PSCOUNT:" + EstimateTransitionMatrix.PSEUDO_COUNTS //
		;
	}

	public EstimateTransitionMatrix(LikelihoodTreeAa tree, MultipleSequenceAlignmentSet msas) {
		this.tree = tree;
		this.msas = msas;
		cacheLogLikelihood = new HashMap<String, Double>();
		numSpecies = tree.childNames().size();
		time = new double[numSpecies][numSpecies];
		N = GprSeq.AMINO_ACIDS.length;
		initNames();
	}

	/**
	 * Calculate 'stable' probability for each amino acid
	 * Note: We calculate using 'all' alignments and the first alignment
	 */
	public ArrayRealVector calcPi() {
		if (pi != null && piVect != null) return piVect;

		System.err.println("Counting amino acids: ");
		int countAa[] = countAa();
		int tot = 0;
		for (int aa = 0; aa < countAa.length; aa++)
			tot += countAa[aa];

		// Label to show results
		String label = "PI_AA";
		if (countAa.length >= 400) label = "PI_AA_PAIR";

		// Calculate for all AA
		double piAll[];
		piAll = new double[countAa.length];
		if (verbose) System.out.println(label + "\tAAcode\tAA\tcount\tp");
		for (int i = 0; i < countAa.length; i++) {
			piAll[i] = countAa[i] / ((double) tot);
			System.out.println(label + "\t" + i + "\t" + names[i] + "\t" + countAa[i] + "\t" + piAll[i]);
		}

		// Create vector
		pi = piAll;
		piVect = new ArrayRealVector(pi);
		return piVect;
	}

	/**
	 * Count occurrences
	 */
	protected int[] countAa() {
		return msas.countAa();
	}

	/**
	 *  Count all transitions
	 */
	protected int[][] countTransitions(int seqNum1, int seqNum2) {
		return msas.countTransitions(seqNum1, seqNum2);
	}

	/**
	 * Inference of a transition matrix Q
	 */
	@SuppressWarnings("unchecked")
	public TransitionMatrix estimateTransitionMatrix() {
		calcPi();

		//----
		// Estimate matrix
		//---

		// For each pair of species, estimate Q
		System.out.println("Estimate transition matrix");
		Counter count = new Counter();
		TransitionMatrix zero = new TransitionMatrix(N, N);

		TransitionMatrix QhatSum = IntStream.range(0, msas.getNumAligns()) //
				.mapToObj(i -> IntStream.range(i + 1, msas.getNumAligns()).mapToObj(j -> new Tuple<Integer, Integer>(i, j))) // Create [i,j] pairs
				.flatMap(s -> s) //
				.map(t -> (Tuple<Integer, Integer>) t) // Cast Object to Tuple
				.parallel() //
				.map(t -> estimateTransitionMatrix(t.first, t.second)) // Estimate transition matrix (can be zero matrix)
				.filter(m -> !m.isZero()) // Remove zero matrices
				.peek(t -> count.inc()) // Count
				.reduce(zero, (a, b) -> a.add(b)) // Add results
		;

		// Calculate the average of all estimators
		Q = new TransitionMatrixMarkov(QhatSum.scalarMultiply(1.0 / count.getCount()));
		Q.setColNames(names);
		Q.setRowNames(names);

		return Q;
	}

	/**
	 * Calculate transition matrix from data
	 */
	public TransitionMatrix estimateTransitionMatrix(int seqNum1, int seqNum2) {
		String seqName1 = msas.getSpecies()[seqNum1];
		String seqName2 = msas.getSpecies()[seqNum2];
		double t = time(seqNum1, seqNum2);
		if (debug) System.err.println(seqName1 + "\t" + seqName2 + "\ttime: " + t);

		//---
		// Count all transitions
		//---
		int count[][] = countTransitions(seqNum1, seqNum2);

		// Add pseudo-counts
		if (PSEUDO_COUNTS > 0) {
			for (int i = 0; i < count.length; i++)
				for (int j = 0; j < count.length; j++)
					count[i][j] += pseudoCount;
		}

		// Calculate total counts
		int sum = Arrays.stream(count).flatMapToInt(x -> Arrays.stream(x)).sum();

		//---
		// Estimate matrix P: Convert to transition frequencies
		//---
		double phat[][] = new double[N][N];
		double n = sum;
		for (int i = 0; i < phat.length; i++) {
			if (pi[i] != 0) {
				long sumRow = 0;
				for (int j = 0; j < phat.length; j++)
					sumRow += count[i][j];

				for (int j = 0; j < phat.length; j++) {
					switch (METHOD) {
					case 0:
						phat[i][j] = count[i][j] / ((double) sumRow);
						break;

					case 1:
						double freq = count[i][j] / n;
						phat[i][j] = freq / pi[i];
						break;

					case 2:
						// We normally use this one
						freq = (count[i][j] + count[j][i]) / (2.0 * n);
						phat[i][j] = freq / pi[i];
						break;

					case 3:
						freq = (count[i][j] + count[j][i]) / (2.0 * n);
						double piAvg = (pi[i] + pi[j]) / 2.0;
						phat[i][j] = freq / piAvg;
						break;

					default:
						throw new RuntimeException("Unimplemented method");
					}
				}
			} else Gpr.debug("WARNING: pi[" + i + "] is zero!");
		}

		// Create transition matrix
		// 		P(t) = exp(t * Q) = V^T exp(t * D) V  => Q = 1/t log[ P(t) ]
		TransitionMatrixMarkov Phat = new TransitionMatrixMarkov(phat);
		TransitionMatrix Qhat = new TransitionMatrixMarkov(Phat.log().scalarMultiply(1 / t));

		// Some sanity checks
		if (Phat.isSymmetric()) Gpr.debug("Phat[" + seqName1 + " , " + seqName2 + "] is symmetric.");
		if (!Phat.isProbabilityMatrix()) Gpr.debug("Phat[" + seqName1 + " , " + seqName2 + "] is NOT a probability matrix.");
		if (Qhat.isSymmetric()) Gpr.debug("Qhat[" + seqName1 + " , " + seqName2 + "] is symmetric.");

		// Remove negative entries from matrix
		if (REMOVE_NEGATIVES > 0) {
			double dqhat[][] = Qhat.getData();
			for (int i = 0; i < dqhat.length; i++)
				for (int j = 0; j < dqhat.length; j++) {
					if (Double.isInfinite(dqhat[i][j]) || Double.isNaN(dqhat[i][j])) throw new RuntimeException("Matrix Qhat contains either NaN or Infinite values: " + seqName1 + ", " + seqName2);
					if (i != j && dqhat[i][j] < 0) dqhat[i][j] = 0;
				}

			// Create matrix
			Qhat = new TransitionMatrixMarkov(dqhat);
		}
		Qhat.setColNames(names);
		Qhat.setRowNames(names);

		// Check
		RealVector z = Qhat.operate(calcPi());
		if (verbose) System.err.println("NORM_QHAT_PI_" + methods() + "\t" + seqName1 + "\t" + seqName2 + "\t" + t + "\t" + z.getNorm());

		return Qhat;
	}

	public RealVector getPi() {
		return piVect;
	}

	public TransitionMatrix getQ() {
		return Q;
	}

	protected void initNames() {
		// Column and row names
		names = new String[GprSeq.AMINO_ACIDS.length];
		for (int i = 0; i < names.length; i++)
			names[i] = GprSeq.code2aa((byte) i) + "";
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
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
