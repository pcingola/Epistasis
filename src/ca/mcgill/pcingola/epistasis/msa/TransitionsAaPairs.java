package ca.mcgill.pcingola.epistasis.msa;

import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.pcingola.epistasis.pdb.DistanceResult;

/**
 * Count amino acid pairs transitions (i.e. 400x400 matrix of AA pairs)
 *
 * @author pcingola
 */
public class TransitionsAaPairs extends TransitionsAa {

	public TransitionsAaPairs() {
		this(true);
	}

	public TransitionsAaPairs(boolean zero) {
		int n = GprSeq.AMINO_ACIDS.length;
		int nn = n * n;
		count = new long[nn][nn];

		seq = new String[nn];
		for (int i = 0, l = 0; i < n; i++)
			for (int j = 0; j < n; j++, l++)
				seq[l] = GprSeq.code2aa((byte) i) + "_" + GprSeq.code2aa((byte) j);

		if (zero) {
			n = count.length;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					count[i][j] = 0;
		}
	}

	public synchronized TransitionsAaPairs add(TransitionsAaPairs t) {
		TransitionsAaPairs sum = new TransitionsAaPairs(false);

		int n = count.length;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				sum.count[i][j] = count[i][j] + t.count[i][j];

		return sum;
	}

	/**
	 * Count sequences transitions (AA pairs)
	 */
	public void count(byte[] codes1, byte[] codes2) {
		int n = GprSeq.AMINO_ACIDS.length;

		for (int i = 0; i < codes1.length; i++) {
			byte c1i = codes1[i], c2i = codes2[i];
			if (c1i < 0 || c2i < 0) continue;
			int idxi = GprSeq.aaPairCode(c1i, c2i);

			for (int j = i + 1; j < codes1.length; j++) {
				byte c1j = codes1[j], c2j = codes2[j];
				if (c1j < 0 || c2j < 0) continue;
				int idxj = GprSeq.aaPairCode(c1j, c2j);

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
}
