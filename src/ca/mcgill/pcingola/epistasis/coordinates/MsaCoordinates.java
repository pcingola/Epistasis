package ca.mcgill.pcingola.epistasis.coordinates;

public class MsaCoordinates {

	public String msaId;

	public int msaIdx;
	public String aaSeq; // Sequence

	public MsaCoordinates(String msaId, int msaIdx) {
		this.msaId = msaId;
		this.msaIdx = msaIdx;
	}

	@Override
	public String toString() {
		return msaId + ":" + msaIdx;
	}
}
