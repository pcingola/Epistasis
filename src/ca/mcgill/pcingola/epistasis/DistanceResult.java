package ca.mcgill.pcingola.epistasis;

import org.biojava.bio.structure.AminoAcid;

public class DistanceResult {

	public String pdbId;
	public String chainId;
	public int aaPos1, aaPos2;
	public char aa1, aa2;
	public double distance;

	DistanceResult(AminoAcid aa1, AminoAcid aa2, double distance) {
		setAa1(aa1);
		setAa2(aa2);
		this.distance = distance;
	}

	public void setAa1(AminoAcid aa) {
		pdbId = aa.getChain().getParent().getPDBCode();
		chainId = aa.getChainId();
		aaPos1 = aa.getResidueNumber().getSeqNum() - 1;
		aa1 = aa.getChemComp().getOne_letter_code().charAt(0);
	}

	public void setAa2(AminoAcid aa) {
		pdbId = aa.getChain().getParent().getPDBCode();
		chainId = aa.getChainId();
		aaPos2 = aa.getResidueNumber().getSeqNum() - 1;
		aa2 = aa.getChemComp().getOne_letter_code().charAt(0);
	}

	@Override
	public String toString() {
		return "pdbId1: " + pdbId //
				+ "\tchain: " + chainId //
				+ "\tdistance: " + distance //
				+ "\taa1: " + aa1 + ", " + aaPos1 //
				+ "\taa2: " + aa2 + ", " + aaPos2 //
		;
	}
}
