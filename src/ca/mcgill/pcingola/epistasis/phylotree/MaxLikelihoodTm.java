package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
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

	boolean verbose = false;

	LikelihoodTree tree;
	MultipleSequenceAlignmentSet msas;
	TransitionMatrixMarkov Q;
	ArrayRealVector piVect;
	double pi[];
	Random random;

	public MaxLikelihoodTm(LikelihoodTree tree, MultipleSequenceAlignmentSet msas) {
		this.tree = tree;
		this.msas = msas;
		random = new Random(20140426);
	}

	/**
	 * Calculate 'stable' probability for each amino acid
	 * Note: We calculate using 'all' alignments and the first alignment
	 */
	protected ArrayRealVector calcPi() {
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
		double t = tree.distance(seqName1, seqName2);
		System.out.println("\t" + seqName1 + "\t" + seqName2 + "\ttime: " + t);

		int pseudoCounts = 1;

		// Count all transitions
		int count[][] = msas.countTransitions(seqNum1, seqNum2);

		// Add pseudo-counts
		for (int i = 0; i < count.length; i++)
			for (int j = 0; j < count.length; j++)
				count[i][j] += pseudoCounts;

		// Calculate total counts
		int sum = Arrays.stream(count).flatMapToInt(x -> Arrays.stream(x)).sum();

		// Convert to transition frequencies
		// Estimate matrix P
		double phat[][] = new double[GprSeq.AMINO_ACIDS.length][GprSeq.AMINO_ACIDS.length];
		double n = sum;
		for (int i = 0; i < phat.length; i++)
			for (int j = 0; j < phat.length; j++)
				phat[i][j] = (count[i][j] + count[j][i]) / n * pi[i];

		// Create matrix
		// P(t) = exp(t * Q) = V^T exp(t * D) V  => Q = 1/t log[ P(t) ]
		TransitionMatrixMarkov Phat = new TransitionMatrixMarkov(phat);
		TransitionMatrix Qhat = new TransitionMatrixMarkov(Phat.log().scalarMultiply(1 / t));
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
	 * Calculate the likelihood of the current transition matrix
	 */
	public double logLikelyhood() {
		double logLik = 0.0;

		// Calculate likelihood for each MSA & each base
		for (MultipleSequenceAlignment msa : msas) {
			System.out.println("MSA: " + msa.getId() + "\t" + logLik);

			for (int i = 0; i < msa.length(); i++) {
				if (msa.isSkip(i)) continue;
				String seqCol = msa.getColumn(i);

				// Set sequence
				tree.setLeafSequence(seqCol);
				double like = tree.likelihood(Q, pi);
				logLik += -Math.log(like);
				if (verbose) System.out.println("Likelyhood: " + like + "\t\t" + logLik + "\t\t" + msa.getId() + "\tpos: " + i + "\t" + seqCol);
			}

		}

		return logLik;
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
			r[i][i] = 50 * random.nextDouble();

		// Only one side
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
				r[i][j] = 0.01 * random.nextDouble();

		// Apply symmetry
		for (int i = 0; i < N; i++)
			for (int j = 0; j < i; j++)
				r[i][j] = r[j][i];

		// Normalize
		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++)
				sum += r[i][j];

			for (int j = 0; j < N; j++)
				r[i][j] /= sum;
		}

		// We assume that the matrix is P(1) = exp( 1.0 * Q )
		TransitionMatrixMarkov P = new TransitionMatrixMarkov(r);
		Gpr.debug("P(1) :\n" + P);

		//---
		// Matrix log
		///---

		// Did we already perform eigendecomposition?
		EigenDecomposition eigen = new EigenDecomposition(P);

		// Exponentiate the diagonal
		RealMatrix D = eigen.getD().copy();
		double maxLambda = Double.NEGATIVE_INFINITY;
		int dim = D.getColumnDimension();
		for (int i = 0; i < dim; i++) {
			double lambda = D.getEntry(i, i);
			maxLambda = Math.max(maxLambda, lambda);
			D.setEntry(i, i, Math.log(lambda));
		}

		// Perform matrix exponential
		RealMatrix qmatrix = eigen.getV().multiply(D).multiply(eigen.getVT());

		Q = new TransitionMatrixMarkov(qmatrix.getData());
		return Q;
	}
}
