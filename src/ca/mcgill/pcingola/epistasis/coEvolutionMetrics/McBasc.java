package ca.mcgill.pcingola.epistasis.coEvolutionMetrics;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Tuple;

/**
 * McBASC: Basic coEvolution algorithm based on correlation between pairs of columns in
 * a multiple sequence alignment
 *
 * References:
 * 		"Correlated mutations and residue contacts in proteins", Goble et. al., 1994
 * 		"Correlated Mutations Contain Information About Protein±protein Interaction", Pazos et.al., 1997
 * 		"Influence of Conservation on Calculations of Amino Acid Covariance in Multiple Sequence Alignments", Fodor et. al., 2004
 *
 * This code is based on Fodor's et. al. implementation, which is based on the other
 * two papers. The similarity matrix is from:
 * 		"Test for comparing related amino acid sequences. Cytochrome C and Cytochrome C551", McLachlan, 1971
 *
 *
 * @author pcingola
 */
public class McBasc {

	/**
	 * When this flag is active, we try to follow Fodor's et. al. implementation.
	 * In their implementation there are few 'particularities'. The authors seem
	 * to be well aware of some of them since they are documented in their code
	 * as bugs or just oddities when trying to follow the specification from the
	 * original papers.
	 */
	protected boolean useFodor = true;

	protected boolean debug = false;
	protected String coli, colj; // Columns 'i' and 'j' from the alignment
	protected byte codei[], codej[]; // Columns transformed from string to codes
	double meanS_i, sigmaS_i;
	double meanS_j, sigmaS_j;
	AaSimilarityMatrix similarytyMatrix;

	/**
	 * Correlation between two columns in an MSA
	 */
	public static double correlation(AaSimilarityMatrix similarytyMatrix, String coli, String colj) {
		McBasc mcBasc = new McBasc(similarytyMatrix, coli, colj);
		mcBasc.setUseFodor(false);
		return mcBasc.score();
	}

	/**
	 * Correlation between two columns in an MSA
	 * Try to follow Fodor's implementation
	 */
	public static double correlationFodor(AaSimilarityMatrix similarytyMatrix, String coli, String colj) {
		McBasc mcBasc = new McBasc(similarytyMatrix, coli, colj);
		mcBasc.setUseFodor(true);
		return mcBasc.score();
	}

	public McBasc(AaSimilarityMatrix similarytyMatrix, String coli, String colj) {
		// Sanity check
		if (coli.length() != colj.length()) throw new RuntimeException("Lengths do not match!");

		// Convert string to byte codes
		this.similarytyMatrix = similarytyMatrix;
		this.coli = coli;
		this.colj = colj;

		int numAligns = coli.length();
		codei = new byte[numAligns];
		codej = new byte[numAligns];

		// Encode columns
		for (int i = 0; i < numAligns; i++) {
			codei[i] = GprSeq.aa2Code(coli.charAt(i));
			codej[i] = GprSeq.aa2Code(colj.charAt(i));
		}
	}

	/**
	 * Calculate <s_i> and sigma_i for a given column
	 */
	Tuple<Double, Double> calcSimilarityStats(byte codei[]) {
		int len = codei.length;

		int count = 0;
		double sum = 0, sum2 = 0;
		for (int k = 0; k < len; k++)
			for (int l = 0; l < len; l++) {

				// Note: Fodor's implementation does use the diagonal to calculating <s_i>, <s_j>,
				//       sigma_i and sigma_j, but does NOT include the diagonal when calculating
				//       the main sum.
				if (((k != l) || useFodor) // Skip diagonal entries?
						&& (codei[k] >= 0) // Skip gaps
						&& (codei[l] >= 0) //
				) {
					if ((codei[k] >= 0) && (codei[l] >= 0)) {
						double sim = similarytyMatrix.getEntry(codei[k], codei[l]);
						sum += sim;
						sum2 += (sim * sim);
						count++;

						if (debug) Gpr.debug("Similarity[ " //
								+ "'" + GprSeq.code2aa(codei[k]) + "': " + codei[k] //
								+ "  , " //
								+ " '" + GprSeq.code2aa(codei[l]) + "': " + codei[l] + " ]" //
								+ ": " + sim //
								+ "\tsum: " + sum //
								+ "\tsum2: " + sum2 //
								+ "\tcount: " + count //
						);
					}
				}
			}
		if (count <= 1) return null;

		double mu = sum / count;

		double n1 = count - 1.0;
		double s2m = (sum * sum) / count;
		double sigma = Math.sqrt((sum2 - s2m) / n1);

		if (debug) Gpr.debug("mu: " + mu + "\tsigma: " + sigma);

		return new Tuple<>(mu, sigma);
	}

	public double getMeanS_i() {
		return meanS_i;
	}

	public double getMeanS_j() {
		return meanS_j;
	}

	public double getSigmaS_i() {
		return sigmaS_i;
	}

	public double getSigmaS_j() {
		return sigmaS_j;
	}

	/**
	 * Calculate McBASC score
	 */
	public double score() {
		// Calculate mean and sigma for column 'i'
		Tuple<Double, Double> msi = calcSimilarityStats(codei);
		if (msi == null) return -1;
		meanS_i = msi.first;
		sigmaS_i = msi.second;
		if (sigmaS_i == 0) return -1; // Perfectly conserved column (except for gaps);

		// Calculate mean and sigma for column 'j'
		Tuple<Double, Double> msj = calcSimilarityStats(codej);
		if (msj == null) return -1;
		meanS_j = msj.first;
		sigmaS_j = msj.second;
		if (sigmaS_j == 0) return -1; // Perfectly conserved column (except for gaps);

		// Main sum
		int N = codei.length;
		int count = 0;
		double sum = 0;
		for (int k = 0; k < N; k++) {
			// Fodor's implementation only adds elements in the upper diagonal
			// This probably has an effect if the matrix is non-symmetric, but
			// it seems OK for McLachlan's matrix.
			int min = useFodor ? k + 1 : 0;

			for (int l = min; l < N; l++) {
				if ((k != l) // Skip diagonal entries? (Note: Fodor's implementation does skip it)
						&& (codei[k] >= 0) && (codei[l] >= 0) //
						&& (codej[k] >= 0) && (codej[l] >= 0) //
				) {
					double s_ikl = similarytyMatrix.getEntry(codei[k], codei[l]);
					double s_jkl = similarytyMatrix.getEntry(codej[k], codej[l]);
					double topi = (s_ikl - meanS_i);
					double topj = (s_jkl - meanS_j);
					double top = topi * topj;
					sum += top;
					count++;
					if (debug) Gpr.debug("k:" + k //
							+ "\tl:" + l //
							+ "\ts_ikl: '" + GprSeq.code2aa(codei[k]) + "' -> '" + GprSeq.code2aa(codei[l]) + "'"//
							+ "\ts_jkl: '" + GprSeq.code2aa(codej[k]) + "' -> '" + GprSeq.code2aa(codej[l]) + "'"//
							+ "\n\t\ts_ikl:" + s_ikl //
							+ "\ts_jkl: " + s_jkl //
							+ "\ttopi: " + topi //
							+ "\ttopj: " + topj //
							+ "\ttop: " + top //
							+ "\tsum: " + sum //
							+ "\tcount:" + count //
					);
				}
			}
		}

		// Correlation coefficient
		if (useFodor) {
			// Note: Fodor's implementation always divides by N^2 independently
			// of the number of term in the sum.
			// Also the first coefficient is 2.0 since it assumes a symmetric
			// matrix and only sums the upper triangle.
			return 2.0 / (sigmaS_i * sigmaS_j * N * N) * sum;
		}
		return 1.0 / (sigmaS_i * sigmaS_j * count) * sum;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setUseFodor(boolean useFodor) {
		this.useFodor = useFodor;
	}

}
