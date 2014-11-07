package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.pdb.PdbGenomeMsas;

/**
 * Store genotype position information and relation to MSA:aaIdx
 *
 * @author pcingola
 */
public class GenotypePos extends Marker {

	private static final long serialVersionUID = 1L;

	protected String msaId = null; // MSA ID information
	protected int aaIdx = -1; // MSA's amino acid index
	protected String annotataions; // Annotations referring to this entry

	public GenotypePos(Marker parent, int start, int end, String id) {
		super(parent, start, end, false, id);
	}

	public GenotypePos(Marker parent, int pos, String id) {
		super(parent, pos, pos, false, id);
	}

	public GenotypePos(String msaId, int aaIdx) {
		super();
		this.msaId = msaId;
		this.aaIdx = aaIdx;
		id = msaId + "[" + aaIdx + "]";
	}

	public GenotypePos(VcfEntry ve) {
		super(ve.getParent(), ve.getStart(), ve.getEnd(), false, ve.getChromosomeName() + ":" + ve.getStart() + "_" + ve.getRef() + "/" + ve.getAltsStr());
		annotataions = ve.getInfo("EFF");
	}

	public int getAaIdx() {
		return aaIdx;
	}

	public String getAnnotataions() {
		return annotataions;
	}

	public String getMsaId() {
		return msaId;
	}

	public boolean hasMsaInfo() {
		return msaId != null;
	}

	/**
	 * Find MSAid and AaIdx for a genomic position (given as an ID string)
	 */
	public boolean mapGenomic2Msa(PdbGenomeMsas pdbGenomeMsas) {
		return mapGenomic2Msa(pdbGenomeMsas, null);
	}

	/**
	 * Find MSAid and AaIdx for a genomic position (given as an ID string)
	 * Map it to transcript 'trId'
	 */
	public boolean mapGenomic2Msa(PdbGenomeMsas pdbGenomeMsas, String trId) {
		// Already mapped? Nothing to do
		if (msaId != null) return true;

		// Create a marker, find all MSAs that intercept the marker
		MultipleSequenceAlignmentSet msas = pdbGenomeMsas.getMsas();
		Markers res = msas.query(this);

		// We now need to find the AA index for that MSA
		for (Marker m : res) {
			String msaId = m.getId();
			MultipleSequenceAlignment msa = msas.getMsa(msaId);

			// Check trancript ID
			if (trId != null && !msa.getTranscriptId().equals(trId)) continue;

			// Map to AA index
			int aaIdx = pdbGenomeMsas.genomicPos2AaIdx(msaId, start);

			// Is AA index within trancript's AA sequence 
			if ((aaIdx >= 0) || (aaIdx >= msa.getAaSeqLen())) {

				Gpr.debug("ERROR: Index out of range !"//
						+ "\n\tID_J              : " + id //
						+ "\n\tMarker            : " + m.toStr() //
						+ "\n\tmsa.Id            : " + msaId //
						+ "\n\tmsa.aaIdx         : " + aaIdx //
				);
			} else {
				// Found it!
				this.msaId = msaId;
				this.aaIdx = aaIdx;
				return true;
			}
		}

		return false;
	}

	/**
	 * Set genomic coordinates based on MsaId + aaIdx
	 * @return true if success
	 */
	public boolean mapMsa2Genomic(PdbGenomeMsas pdbGenomeMsas) {
		return mapMsa2GenomicErr(pdbGenomeMsas) == null;
	}

	/**
	 * Set genomic coordinates based on MsaId + aaIdx
	 * @return String describing error type, null on success
	 */
	public String mapMsa2GenomicErr(PdbGenomeMsas pdbGenomeMsas) {
		// Already mapped? Nothing to do
		if (start >= 0) return null;

		if (aaIdx < 0) return "AA index is negative"; // Incorrect aaIdx

		// Find MSA
		MultipleSequenceAlignmentSet msas = pdbGenomeMsas.getMsas();
		MultipleSequenceAlignment msa = msas.getMsa(msaId);
		if (msa == null) return "MSA '" + msaId + "' not found";

		// Get transcript
		String trId = msa.getTranscriptId();
		Transcript tr = pdbGenomeMsas.getTranscript(trId);
		if (tr == null) return "Transcript '" + trId + "' not found";

		// Find genomic position based on AA position
		int aa2pos[] = tr.aaNumber2Pos();
		if (aa2pos.length <= aaIdx) return "AA index out of range (aaIdx: " + aaIdx + ", protein length: " + aa2pos.length + ")"; // aaIdx Outside transcript

		// Convert to genomic positions
		start = end = aa2pos[aaIdx];
		parent = tr.getChromosome();

		return null;
	}

	public void setMsa(String msaId, int aaIDx) {
		this.msaId = msaId;
		aaIdx = aaIDx;
	}

}
