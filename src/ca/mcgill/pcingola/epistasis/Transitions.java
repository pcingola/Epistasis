package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Count transision pairs (i.e. 400x400 matrix of AA pairs)
 *
 * @author pcingola
 */
public class Transitions {

	String seq[];
	long count[][];

	public Transitions() {
		int n = GprSeq.AMINO_ACIDS.length;
		int nn = n * n;
		count = new long[nn][nn];

		seq = new String[nn];
		for (int i = 0, l = 0; i < n; i++)
			for (int j = 0; j < n; j++, l++)
				seq[l] = GprSeq.code2aa((byte) i) + "_" + GprSeq.code2aa((byte) j);

	}

	/**
	 * Count sequences transitions (AA pairs)
	 */
	public void count(byte[] codes1, byte[] codes2) {
		int n = GprSeq.AMINO_ACIDS.length;

		for (int i = 0; i < codes1.length; i++) {
			byte c1i = codes1[i], c2i = codes2[i];
			if (c1i < 0 || c2i < 0) continue;
			int idxi = c1i * n + c2i;

			for (int j = i + 1; j < codes1.length; j++) {
				byte c1j = codes1[j], c2j = codes2[j];
				if (c1j < 0 || c2j < 0) continue;
				int idxj = c1j * n + c2j;

				// Symmetric count
				count[idxi][idxj]++;
				count[idxj][idxi]++;
			}
		}
	}

	/**
	 * Count sequences transitions (AA pairs)
	 */
	public void count(DistanceResult d) {
		count(GprSeq.aa2Code(d.aaSeq1), GprSeq.aa2Code(d.aaSeq2));
	}

	public long[][] getCount() {
		return count;
	}

	/**
	 * Calculate 'row normalized' frequencies
	 * @return
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
		long sum = 0;
		for (int i = 0; i < count.length; i++) {
			sb.append(seq[i]);
			for (int j = 0; j < count.length; j++) {
				sb.append("\t" + count[i][j]);
				sum += count[i][j];
			}
			sb.append("\n");
		}

		sb.append("Sum: " + sum);
		return sb.toString();
	}

	public String toString(double f[][]) {
		StringBuilder sb = new StringBuilder();

		// Title
		for (int i = 0; i < f.length; i++)
			sb.append("\t" + seq[i]);
		sb.append("\n");

		// Values
		double sum = 0;
		for (int i = 0; i < f.length; i++) {
			sb.append(seq[i]);
			for (int j = 0; j < f.length; j++) {
				if (f[i][j] == 0.0) sb.append("\t0");
				else if (f[i][j] == 1.0) sb.append("\t1");
				else sb.append(String.format("\t%.3e", f[i][j]));
				sum += f[i][j];
			}
			sb.append("\n");
		}

		sb.append("Sum: " + sum);
		return sb.toString();
	}
}
