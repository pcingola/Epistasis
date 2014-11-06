package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
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
	public boolean map2MsaAa(PdbGenomeMsas pdbGenomeMsas) {

		// Already mapped? Nothing to do
		if (msaId != null) return true;

		// Create a marker, find all MSAs that intercept the marker
		MultipleSequenceAlignmentSet msas = pdbGenomeMsas.getMsas();
		Markers res = msas.query(this);

		// We now need to find the AA index for that MSA
		for (Marker m : res) {
			String msaId = m.getId();
			int aaIdx = pdbGenomeMsas.genomicPos2AaIdx(msaId, start);

			if (aaIdx >= 0) {
				// Found something: Store result
				// Note: We could just stop here, all the rest is done just to
				//       ensure that the mapping works OK and that we are not
				//       finding inconsistent sequences

				// Check sequence length
				MultipleSequenceAlignment msa = msas.getMsa(msaId);
				if (aaIdx >= msa.getSeqLen()) {
					Gpr.debug("ERROR: Index out of range !"//
							+ "\n\tID_J              : " + id //
							+ "\n\tMarker            : " + m.toStr() //
							+ "\n\tmsa.Id            : " + msaId //
							+ "\n\tmsa.aaIdx         : " + aaIdx //
							);
				} else {
					this.msaId = msaId;
					this.aaIdx = aaIdx;
					return true;
				}
			}
		}

		return false;
	}

	public void setMsa(String msaId, int aaIDx) {
		this.msaId = msaId;
		aaIdx = aaIDx;
	}

}
