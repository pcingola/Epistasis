package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignment;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;
import ca.mcgill.pcingola.epistasis.pdb.PdbGenomeMsas;

/**
 * Store genotype information
 *
 * @author pcingola
 */
public class Genotype extends Marker {

	private static final long serialVersionUID = 1L;

	int minorAlleleCount;
	byte gt[];
	String msaId = null; // MSA ID information
	int aaIdx = -1; // MSA's amino acid index
	String annotataions; // Annotations referring to this entry

	public Genotype(Marker parent, int start, int end, String id, byte gt[]) {
		super(parent, start, end, false, id);
		this.gt = minorAllele(gt);
	}

	public Genotype(VcfEntry ve) {
		super(ve.getParent(), ve.getStart(), ve.getEnd(), false, ve.getChromosomeName() + ":" + ve.getStart() + "_" + ve.getRef() + "/" + ve.getAltsStr());
		gt = minorAllele(ve.getGenotypesScores());
		annotataions = ve.getInfo("EFF");
	}

	public int getAaIdx() {
		return aaIdx;
	}

	public byte[] getGt() {
		return gt;
	}

	public int getMinorAlleleCount() {
		return minorAlleleCount;
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
	public void map2MsaAa(PdbGenomeMsas pdbGenomeMsas) {

		// Already mapped? Nothing to do
		if (msaId != null) return;

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
					return;
				}
			}
		}
	}

	/**
	 * Convert to minor allele (or filter out)
	 * @return A minor allele genotype, or null if it doesn't satisfy some fitlering requirements)
	 */
	byte[] minorAllele(byte gt[]) {
		// Count alleles
		minorAlleleCount = 0;
		for (int i = 0; i < gt.length; i++)
			if (gt[i] > 0) minorAlleleCount += gt[i]; // Don't count '-1' (i.e. missing genotypes)

		if (minorAlleleCount <= gt.length) return gt; // OK, gt[] is mainor allele

		// Convert to minor allele
		for (int i = 0; i < gt.length; i++)
			if (gt[i] >= 0) gt[i] = (byte) (2 - gt[i]);

		// Convert to minor allele
		minorAlleleCount = 2 * gt.length - minorAlleleCount;

		return gt;

	}

	public void setMsa(String msaId, int aaIDx) {
		this.msaId = msaId;
		aaIdx = aaIDx;
	}

}
