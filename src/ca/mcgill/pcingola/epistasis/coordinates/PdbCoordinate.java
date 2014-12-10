package ca.mcgill.pcingola.epistasis.coordinates;

public class PdbCoordinate {

	public String pdbId;
	public String pdbChainId;
	public int aaPos;
	public char aa;

	@Override
	public String toString() {
		return pdbId + ":" + pdbChainId + "[" + aaPos + "] '" + aa + "'";
	}
}
