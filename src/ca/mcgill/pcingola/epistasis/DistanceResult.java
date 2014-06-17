package ca.mcgill.pcingola.epistasis;

import org.biojava.bio.structure.AminoAcid;

import ca.mcgill.mcb.pcingola.util.Gpr;

public class DistanceResult {

	public String pdbId;
	public String pdbChainId;
	public int aaPos1, aaPos2;
	public char aa1, aa2;
	public double distance;
	public String aaSeq1, aaSeq2;
	public String chr1, chr2;
	public int pos1, pos2;

	DistanceResult(AminoAcid aa1, AminoAcid aa2, double distance) {
		setAa1(aa1);
		setAa2(aa2);
		this.distance = distance;
	}

	DistanceResult(String line) {
		String fields[] = line.split("\t");
		int n = 0;
		pdbId = fields[n++];
		pdbChainId = fields[n++];
		distance = Gpr.parseDoubleSafe(fields[n++]);
		aa1 = fields[n++].charAt(0);
		aaPos1 = Gpr.parseIntSafe(fields[n++]);
		aa2 = fields[n++].charAt(0);
		aaPos2 = Gpr.parseIntSafe(fields[n++]);

		// Optional fields
		if (fields.length > n) {
			String chrPos1 = fields[n++];
			String f[] = chrPos1.split(":");
			chr1 = f[0];
			pos1 = Gpr.parseIntSafe(f[1]);

			String chrPos2 = fields[n++];
			f = chrPos2.split(":");
			chr2 = f[0];
			pos2 = Gpr.parseIntSafe(f[1]);

			aaSeq1 = fields[n++];
			aaSeq2 = fields[n++];
		}
	}

	/**
	 * Compare by genomic position
	 */
	public int compareByPos(DistanceResult d) {
		// Compare first position
		int comp = chr1.compareTo(d.chr1);
		if (comp != 0) return comp;

		comp = pos1 - d.pos1;
		if (comp != 0) return comp;

		// Compare second position
		comp = chr2.compareTo(d.chr2);
		if (comp != 0) return comp;

		comp = pos2 - d.pos2;
		if (comp != 0) return comp;

		// Compare distances
		return (int) Math.signum(distance - d.distance);
	}

	/**
	 * Same genomic positions
	 */
	public boolean equalPos(DistanceResult d) {
		return chr1.equals(d.chr1) //
				&& chr2.equals(d.chr2) //
				&& pos1 == d.pos1 //
				&& pos2 == d.pos2 //
				;
	}

	public void setAa1(AminoAcid aa) {
		pdbId = aa.getChain().getParent().getPDBCode();
		pdbChainId = aa.getChainId();
		aaPos1 = aa.getResidueNumber().getSeqNum() - 1;
		aa1 = aa.getChemComp().getOne_letter_code().charAt(0);
	}

	public void setAa2(AminoAcid aa) {
		pdbId = aa.getChain().getParent().getPDBCode();
		pdbChainId = aa.getChainId();
		aaPos2 = aa.getResidueNumber().getSeqNum() - 1;
		aa2 = aa.getChemComp().getOne_letter_code().charAt(0);
	}

	@Override
	public String toString() {
		return pdbId //
				+ "\t" + pdbChainId //
				+ "\t" + distance //
				+ "\t" + aa1 //
				+ "\t" + aaPos1 //
				+ "\t" + aa2 //
				+ "\t" + aaPos2 //
				+ (chr1 != null ? "\t" + chr1 + ":" + pos1 : "") //
				+ (chr2 != null ? "\t" + chr2 + ":" + pos2 : "") //
				+ (aaSeq1 != null ? "\t" + aaSeq1 : "") //
				+ (aaSeq2 != null ? "\t" + aaSeq2 : "") //
				;
	}

	public String toStringPos() {
		return "" //
				+ (chr1 != null ? "\t" + chr1 + ":" + pos1 : "") //
				+ (chr2 != null ? "\t" + chr2 + ":" + pos2 : "") //
				;
	}
}
