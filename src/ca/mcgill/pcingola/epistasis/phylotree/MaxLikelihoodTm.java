package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

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
		int countAaFirst[] = msas.countAa(0); // Count AA only using the first sequence
		int tot = 0, totFirst = 0;
		for (int aa = 0; aa < countAa.length; aa++) {
			tot += countAa[aa];
			totFirst += countAaFirst[aa];
		}

		// Calculate for all AA
		double piAll[], piFirst[];
		piAll = new double[countAa.length];
		piFirst = new double[countAa.length];
		if (verbose) System.out.println("AAcode\tAA\tcount_all\tcount_first\tp_all\tp_first");
		for (int aa = 0; aa < countAa.length; aa++) {
			piAll[aa] = countAa[aa] / ((double) tot);
			piFirst[aa] = countAaFirst[aa] / ((double) totFirst);
			if (verbose) System.out.println(aa + "\t" + GprSeq.code2aa((byte) aa) + "\t" + countAa[aa] + "\t" + countAaFirst[aa] + "\t" + piAll[aa] + "\t" + piFirst[aa]);
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

		// For each pair of species, estimate Q
		System.out.println("Estimate transition matrix");
		int N = GprSeq.AMINO_ACIDS.length;
		Array2DRowRealMatrix QhatSum = new Array2DRowRealMatrix(N, N);
		int count = 0;
		for (int i = 0; i < msas.getNumAligns(); i++) {
			for (int j = i + 1; j < msas.getNumAligns(); j++) {
				TransitionMatrix QhatTmp = estimateTransitionMatrix(i, j);
				QhatSum = QhatSum.add(QhatTmp);
				count++;
			}
		}

		// Calculate the average of all estimators
		Q = new TransitionMatrixMarkov(QhatSum.scalarMultiply(1.0 / count));
		System.out.println("Qhat: " + count + " estimations\n" + Q);

		return Q;
	}

	/**
	 * Calculate transition matrix from data
	 * @param seqNum1
	 * @param seqNum2
	 */
	public TransitionMatrix estimateTransitionMatrix(int seqNum1, int seqNum2) {
		String seqName1 = msas.getSpecies()[seqNum1];
		String seqName2 = msas.getSpecies()[seqNum2];
		double t = time(seqNum1, seqNum2);
		System.out.println("\t" + seqName1 + "\t" + seqName2 + "\ttime: " + t);

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
			for (int j = 0; j < phat.length; j++)
				phat[i][j] = (count[i][j] + count[j][i]) / (n * pi[i]); // Note: We use symmetry

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
	 * @param msa
	 * @param pos
	 * @return
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
	 * Create a random transition matrix
	 *
	 * Note: We make sure the matrix satisfies symmetry, detailed balance
	 * and stability (the latter is implied by Gershgorin's theorem from
	 * Markov's detailed balance condition)
	 */
	public TransitionMatrix randQ() {
		int N = GprSeq.AMINO_ACIDS.length;
		double r[][] = new double[N][N];

		if (pi == null) calcPi();

		//---
		// Random symmetric probability matrix
		//---

		// Diagonal
		for (int i = 0; i < N; i++)
			r[i][i] = 0 * random.nextDouble();

		// Only one side
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
				r[i][j] = random.nextDouble();

		// Apply symmetry
		for (int i = 0; i < N; i++)
			for (int j = 0; j < i; j++)
				r[i][j] = r[j][i];

		// Normalize
		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++)
				sum += r[i][j];

			r[i][i] = -sum;
		}

		// We assume that the matrix is P(1) = exp( 1.0 * Q )
		// TransitionMatrixMarkov P = new TransitionMatrixMarkov(r);
		// Gpr.debug("P (1):\n" + P);
		Q = new TransitionMatrixMarkov(r);

		//---
		// Matrix log
		///---
		//		// Did we already perform eigendecomposition?
		//		EigenDecomposition eigen = new EigenDecomposition(P);
		//
		//		// Exponentiate the diagonal
		//		RealMatrix D = eigen.getD().copy();
		//		double maxLambda = Double.NEGATIVE_INFINITY;
		//		int dim = D.getColumnDimension();
		//		for (int i = 0; i < dim; i++) {
		//			double lambda = D.getEntry(i, i);
		//			maxLambda = Math.max(maxLambda, lambda);
		//			Gpr.debug("\tLambda: " + lambda + "\tMax lambda: " + maxLambda);
		//			// D.setEntry(i, i, Math.log(lambda));
		//		}
		//		Gpr.debug("Max lambda: " + maxLambda);
		//
		//		// Perform matrix exponential
		//		RealMatrix qmatrix = eigen.getV().multiply(D).multiply(eigen.getVT());
		//
		//		Q = new TransitionMatrixMarkov(qmatrix.getData());
		return Q;
	}

	/**
	 * Show matrix's eigenvalues
	 */
	public void showEienQ() {
		// Did we already perform eigendecomposition?
		EigenDecomposition eigen = new EigenDecomposition(Q);

		// Exponentiate the diagonal
		System.out.println("Q's Eigenvalues: ");
		RealMatrix D = eigen.getD().copy();
		double maxLambda = Double.NEGATIVE_INFINITY;
		int dim = D.getColumnDimension();
		for (int i = 0; i < dim; i++) {
			double lambda = D.getEntry(i, i);
			maxLambda = Math.max(maxLambda, lambda);
			System.out.println("\tlambda_" + i + ":\t" + lambda);
		}

		System.out.println("\tlambda_max:\t" + maxLambda);
		if (maxLambda > 0) throw new RuntimeException("All Q's eigenvalues should be non-positive!");
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
