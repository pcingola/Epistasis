package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
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

	public static int METHOD = 1;

	boolean verbose = false;
	boolean debug = false;

	int N;
	int pseudoCount = 1;
	int numSpecies;
	double pi[];
	double time[][];
	String names[];
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
		N = GprSeq.AMINO_ACIDS.length;
		initNames();
	}

	/**
	 * Calculate 'stable' probability for each amino acid
	 * Note: We calculate using 'all' alignments and the first alignment
	 */
	public ArrayRealVector calcPi() {
		if (pi != null && piVect != null) return piVect;

		System.out.println("Counting amino acids: ");
		int countAa[] = countAa();
		int tot = 0;
		for (int aa = 0; aa < countAa.length; aa++)
			tot += countAa[aa];

		// Calculate for all AA
		double piAll[];
		piAll = new double[countAa.length];
		if (verbose) System.out.println("AAcode\tAA\tcount_all\tcount_first\tp_all\tp_first");
		for (int i = 0; i < countAa.length; i++) {
			piAll[i] = countAa[i] / ((double) tot);
			if (verbose) System.out.println(i + "\t" + names[i] + "\t" + countAa[i] + "\t" + piAll[i]);
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
				.mapToObj(i -> IntStream.range(i + 1, msas.getNumAligns()).mapToObj(j -> new Tuple<Integer,Integer>(i,j))) //
				.flatMap(s -> s) //
				.map(t -> (Tuple<Integer,Integer>)t) //
				.parallel() //
				.map( t -> estimateTransitionMatrix(t.first, t.second) ) //
				.peek( t-> count.inc() ) //
				.reduce( zero, (a,b) -> a.add(b) );
		;

//		Array2DRowRealMatrix QhatSum = new Array2DRowRealMatrix(N, N);
//				for (int i = 0; i < msas.getNumAligns(); i++) {
//					for (int j = i + 1; j < msas.getNumAligns(); j++) {
//						TransitionMatrix QhatTmp = estimateTransitionMatrix(i, j);
//		
//						// Add all transition matrix estimates
//						QhatSum = QhatSum.add(QhatTmp);
//						count++;
//					}
//				}

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
		if(debug) System.err.println(seqName1 + "\t" + seqName2 + "\ttime: " + t);

		// Count all transitions
		int count[][] = countTransitions(seqNum1, seqNum2);

		// Add pseudo-counts
		for (int i = 0; i < count.length; i++) {
			long sumRow=0;
			for (int j = 0; j < count.length; j++) {
				sumRow += count[i][j];
				count[i][j] += pseudoCount;
			}
			if( sumRow==0) Gpr.debug("Row sum is zero! Species " + seqName1 + ", " + seqName2 + ", row "+i);
		}

		// Calculate total counts
		int sum = Arrays.stream(count).flatMapToInt(x -> Arrays.stream(x)).sum();

		// Convert to transition frequencies
		// Estimate matrix P
		double phat[][] = new double[N][N];
		double n = sum;
		for (int i = 0; i < phat.length; i++) {
			if (pi[i] != 0) {
				for (int j = 0; j < phat.length; j++) {
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
			}
		}

		// Create matrix
		// P(t) = exp(t * Q) = V^T exp(t * D) V  => Q = 1/t log[ P(t) ]
		TransitionMatrixMarkov Phat = new TransitionMatrixMarkov(phat);
		TransitionMatrix Qhat = new TransitionMatrixMarkov(Phat.log().scalarMultiply(1 / t));

		// Remove negative entries from matrix
		double dqhat[][] = Qhat.getData();
		for (int i = 0; i < dqhat.length; i++)
			for (int j = 0; j < dqhat.length; j++) {
				if (Double.isInfinite(dqhat[i][j]) || Double.isNaN(dqhat[i][j])) {
					Gpr.toFile("Count." + seqName1 + "_" + seqName2 + ".txt", new TransitionMatrix(count));
					Gpr.toFile("Phat." + seqName1 + "_" + seqName2 + ".txt", Phat.toString());
					Gpr.toFile("Qhat." + seqName1 + "_" + seqName2 + ".txt", Qhat.toString());
					throw new RuntimeException("Matrix Qhat contains either NaN or Infinite values: " + seqName1 + ", " + seqName2);
				}
				if (i != j && dqhat[i][j] < 0) dqhat[i][j] = 0;
			}

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

	protected void initNames() {
		// Column and row names
		names = new String[GprSeq.AMINO_ACIDS.length];
		for (int i = 0; i < names.length; i++)
			names[i] = GprSeq.code2aa((byte) i) + "";
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
