package ca.mcgill.pcingola.epistasis;

import org.biojava.bio.structure.AminoAcid;

import ca.mcgill.mcb.pcingola.util.Gpr;

public class DistanceResult {

	public String pdbId;
	public String pdbChainId;
	public int aaPos1, aaPos2;
	public char aa1, aa2;
	public double distance;

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
				;
	}
}
