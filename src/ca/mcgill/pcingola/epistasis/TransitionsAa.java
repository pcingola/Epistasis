package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Count amino acid transition (i.e. 20x20 matrix of AA)
 *
 * @author pcingola
 */
public class TransitionsAa {

	String seq[];
	long count[][];

	public TransitionsAa() {
		this(true);
	}

	public TransitionsAa(boolean zero) {
		int n = GprSeq.AMINO_ACIDS.length;
		count = new long[n][n];

		seq = new String[n];
		for (int i = 0; i < n; i++)
			seq[i] = GprSeq.code2aa((byte) i) + "";

		if (zero) {
			n = count.length;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					count[i][j] = 0;
		}
	}

	/**
	 * Add a transition matrix
	 */
	public synchronized TransitionsAa add(TransitionsAa t) {
		TransitionsAa sum = new TransitionsAa(false);

		int n = count.length;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				sum.count[i][j] = count[i][j] + t.count[i][j];

		return sum;
	}

	/**
	 * Count transitions
	 */
	public void count(byte[] codes) {
		int n = codes.length;

		for (int i = 0; i < n; i++) {
			int idxi = codes[i];
			if (idxi < 0) continue;

			for (int j = i + 1; j < n; j++) {
				int idxj = codes[j];
				if (idxj < 0) continue;

				// Symmetric count
				count[idxi][idxj]++;
				count[idxj][idxi]++;
			}
		}
	}

	public long[][] getCount() {
		return count;
	}

	/**
	 * Calculate 'row normalized' frequencies
	 */
	public double[][] getFrequenciesRow() {
		int n = count.length;
		double f[][] = new double[n][n];

		for (int i = 0; i < n; i++) {
			// Row sum
			long sum = 0;
			for (int j = 0; j < n; j++)
				sum += count[i][j];

			// Row frequencies
			for (int j = 0; j < n; j++)
				if (sum > 0) f[i][j] = ((double) count[i][j]) / ((double) sum);
				else f[i][j] = 0;
		}

		return f;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		// Title
		for (int i = 0; i < count.length; i++)
			sb.append((i > 0 ? "\t" : "") + seq[i]);
		sb.append("\n");

		// Values
		for (int i = 0; i < count.length; i++) {
			sb.append(seq[i]);
			for (int j = 0; j < count.length; j++)
				sb.append("\t" + count[i][j]);
			sb.append("\n");
		}
		return sb.toString();
	}

	public String toString(double f[][]) {
		StringBuilder sb = new StringBuilder();

		// Title
		for (int i = 0; i < f.length; i++)
			sb.append("\t" + seq[i]);
		sb.append("\n");

		// Values
		for (int i = 0; i < f.length; i++) {
			sb.append(seq[i]);
			for (int j = 0; j < f.length; j++) {
				if (f[i][j] == 0.0) sb.append("\t0");
				else if (f[i][j] == 1.0) sb.append("\t1");
				else sb.append(String.format("\t%.3e", f[i][j]));
			}
			sb.append("\n");
		}

		return sb.toString();
	}
}
