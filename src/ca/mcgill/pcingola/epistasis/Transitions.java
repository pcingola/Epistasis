package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.GprSeq;

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
				seq[l] = "" + GprSeq.code2aa((byte) i) + GprSeq.code2aa((byte) j);

	}

	/**
	 * Count sequences transitions
	 * @param d
	 */
	public void count(DistanceResult d) {
		int n = GprSeq.AMINO_ACIDS.length;
		byte[] codes1 = GprSeq.aa2Code(d.aaSeq1);
		byte[] codes2 = GprSeq.aa2Code(d.aaSeq2);

		boolean ord = (codes1[0] <= codes2[0]);

		for (int i = 0; i < codes1.length; i++) {
			byte c1i = codes1[i], c2i = codes2[i];
			if (c1i < 0 || c2i < 0) continue;
			int idxi = (ord ? c1i * n + c2i : c2i * n + c1i);

			for (int j = i + 1; j < codes1.length; j++) {
				byte c1j = codes1[j], c2j = codes2[j];
				if (c1j < 0 || c2j < 0) continue;
				int idxj = (ord ? c1j * n + c2j : c2j * n + c1j);

				count[idxi][idxj]++;
				if (idxi != idxj) count[idxj][idxi]++;
			}
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();

		// Title
		for (int i = 0; i < count.length; i++)
			sb.append("\t" + seq[i]);
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
}
