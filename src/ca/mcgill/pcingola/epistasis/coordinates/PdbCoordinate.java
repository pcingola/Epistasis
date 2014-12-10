package ca.mcgill.pcingola.epistasis.coordinates;

public class PdbCoordinate {

	public String pdbId;
	public String pdbChainId;
	public int aaPos;
	public char aa;

	public PdbCoordinate(String pdbId, String pdbChainId, int aaPos) {
		this(pdbId, pdbChainId, aaPos, ' ');
	}

	public PdbCoordinate(String pdbId, String pdbChainId, int aaPos, char aa) {
		this.pdbId = pdbId;
		this.pdbChainId = pdbChainId;
		this.aaPos = aaPos;
		this.aa = aa;
	}

	@Override
	public String toString() {
		return pdbId + ":" + pdbChainId + "[" + aaPos + "]" + (aa > 0 ? " '" + aa + "'" : "");
	}
}
