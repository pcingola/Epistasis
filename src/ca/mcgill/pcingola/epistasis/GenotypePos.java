package ca.mcgill.pcingola.epistasis;

import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Exon;
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

	public static boolean debug = false;

	protected String msaId = null; // MSA ID information
	protected int aaIdx = -1; // MSA's amino acid index. Note that aaIndex is respect to an MSA (which is an alignment of a single exon), and not respect to the whole AA sequence
	protected String annotataions; // Annotations referring to this entry

	public GenotypePos() {
		super();
	}

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
			int aaIdx = mapMsaIdPos2AaIdx(pdbGenomeMsas, msaId, start);

			// Is AA index within trancript's AA sequence
			if ((aaIdx < 0) || (aaIdx >= msa.getAaSeqLen())) {
				Gpr.debug("ERROR: Index out of range !"//
						+ "\n\tID                : " + id //
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
		if (!tr.intersects(msa)) return "Transcript and MSA do not intersect: " + tr.toStr() + ", " + msa.toStr();

		// Find genomic position based on AA position
		int aa2pos[] = tr.aaNumber2Pos();

		// Find an exon that matches MSA
		Exon exon = tr.findExon(msa);
		if (exon == null) return "Could not find exon intersecting MSA coordinates";

		// Index out of exon boundaries?
		int exAaEnd = exon.getAaIdxEnd();
		int exAaStart = exon.getAaIdxStart();
		if (exon.getFrame() == 1) exAaStart--; // I don't know why UCSC numbers the AA different when frame is 2
		int numAasEx = exAaEnd - exAaStart + 1;
		if (aaIdx > numAasEx) return "AA index out of range. AaIdx: " + aaIdx + ", exon's AA length: " + numAasEx;

		// Calculate absolute AA index (in full AA sequence)
		int aaSeqIndex = aaIdx + exAaStart;

		// The codon for this amino acid is split by an intron
		int pos;
		if (aaIdx == numAasEx) pos = aa2pos[aaSeqIndex - 1] + 3;
		else if ((aaIdx == 0) && (exon.getFrame() == 1)) pos = exon.isStrandPlus() ? exon.getStart() : exon.getEnd();
		else pos = aa2pos[aaSeqIndex];

		// Convert to AA sequence index to genomic position
		start = end = pos;
		parent = tr.getChromosome();

		// Does this position match MSA coordinates?
		if (!msa.intersects(this)) return "Calculated genomic positions '" + tr.getChromosomeName() + ":" + start + "' in not included in MSA coordinates " + msa.getChromosomeName() + ":" + msa.getStart() + "-" + msa.getEnd();

		// Check AA sequence
		String protein = tr.protein();
		if (aaSeqIndex > protein.length()) return "AA sequence index (aaSeqIndex: " + aaSeqIndex + ") outside AA sequence (length: " + protein.length() + ")";

		// Different AA sequence?
		char aaTr = protein.charAt(aaSeqIndex);
		char aaMsa = msa.getChar(0, aaIdx);
		if (aaTr != aaMsa) {
			if ((aaTr == '*') && (aaMsa == '-')) {
				// Stop codons are represented using different chars, so this one is OK
			} else {
				// All other ones are considered errors
				return "AA from MSA ('" + msa.getChar(0, aaIdx) + "') does not match AA from protein ('" + protein.charAt(aaSeqIndex) + "')";
			}
		}

		// OK, no errors
		return null;
	}

	/**
	 * Get aaIdx from MSA and genomic position  'pos'
	 * Note: Chromosme name is implied by msaId
	 */
	int mapMsaIdPos2AaIdx(PdbGenomeMsas pdbGenomeMsas, String msaId, int pos) {
		// Find all MSA
		MultipleSequenceAlignment msa = pdbGenomeMsas.getMsas().getMsa(msaId);
		if (msa == null) return -1;

		String trid = msa.getTranscriptId();
		Transcript tr = pdbGenomeMsas.getTranscript(trid);

		// Return column index
		return mapMsaTrPos2AaIdx(msa, tr, pos);
	}

	/**
	 * Map a genomic position to an MSA index (given the Transcript and the MSA)
	 * @return -1 on failure
	 */
	int mapMsaTrPos2AaIdx(MultipleSequenceAlignment msa, Transcript tr, int pos) {
		if (tr == null) return -1;

		// Check all MSA
		// Different chromosome or position? Skip
		if (!msa.intersects(tr)) return -1;

		// Find exon
		Exon exon = tr.findExon(pos);
		if (exon == null) {
			Gpr.debug("Cannot find exon for position " + pos + " in transcript " + tr.getId());
			return -1;
		}

		// Find index
		int idxBase = tr.isStrandPlus() ? (pos - msa.getStart()) : (msa.getEnd() - pos);
		int idxAa = idxBase / 3;

		// WARNIGN: If exon frame is 1, the MSA has one additional AA (from the previous exon).
		//          I don't know why they do it this way...
		if (exon.getFrame() == 1) {
			if (idxBase < 1) idxAa = 0; // First two bases are AA number zero
			else idxAa++; // Other bases are AA number 1 and on
		}

		// Return column index
		return idxAa;
	}

	/**
	 * Convert <transcript_id, position> to <msaIdx, aaIdx>
	 */
	public boolean mapTrPos2MsaIdx(PdbGenomeMsas pdbGenomeMsas, String trid, int pos) {
		// Already mapped?
		if (msaId != null) return true;

		// Find transcript
		Transcript tr = pdbGenomeMsas.getTranscript(trid);
		if (tr == null) return false;

		// Set genomic coordinates
		parent = tr.getChromosome();
		start = end = pos;

		// Find all MSAs for a given transcript ID
		List<MultipleSequenceAlignment> msaList = pdbGenomeMsas.getMsas().getMsasByTrId(trid);
		if (msaList == null) return false;

		// Check all MSA
		for (MultipleSequenceAlignment msa : msaList) {
			// Does this MSA intersect chr:pos?
			if (!msa.intersects(this)) continue;

			// Find aaIndx
			int idxAa = mapMsaTrPos2AaIdx(msa, tr, pos);
			if (idxAa < 0) return false;

			// OK
			aaIdx = idxAa;
			msaId = msa.getId();
			return true;
		}

		return false;
	}

	/**
	 * Set this marker to encompass an amino acid (trId:aaIdx)
	 *
	 * Important: 	The marker has ALL bases in the codon.
	 * 				For instance, is the codon is split between two exons, the
	 * 				marker will contain the intron

	 * @return true if successful
	 */
	public boolean markerTrAaIdx(PdbGenomeMsas pdbGenomeMsas, String trId, int aaIdx, char aaExpected) {
		// Find transcript and exon
		Transcript tr = pdbGenomeMsas.getTranscript(trId);
		if (tr == null) return false;

		Exon ex = tr.findExon(start);
		if (ex == null) return false;

		// Calculate start position
		int startPos;
		int fr = 0;
		if (ex.getFrame() != 0) {
			aaIdx--;
			if (ex.getFrame() == 2) aaIdx++; // I don't know why UCSC numbers the AA different when frame is 2
			fr = 3 - ex.getFrame(); // Offset based on frame
		}

		// Find AA start position
		if (ex.isStrandPlus()) {
			int exStart = Math.max(start, tr.getCdsStart());
			startPos = exStart + (aaIdx * 3 + fr);
		} else {
			int exEnd = Math.min(end, tr.getCdsStart());
			startPos = exEnd - (aaIdx * 3 + fr);
		}

		// Get position within CDS
		int cdsBase = tr.baseNumberCds(startPos, false);
		int cds2pos[] = tr.baseNumberCds2Pos();
		if ((ex.isStrandPlus() && (startPos < ex.getStart())) //
				|| (ex.isStrandMinus() && (startPos > ex.getEnd()))) {
			// If the position is outside the exon, then we must jump to previous exon
			startPos = cds2pos[cdsBase - ex.getFrame()];
			cdsBase = tr.baseNumberCds(startPos, true);
		}

		//---
		// Sanity check: Make sure that AA matches between transcript model and MSA data from 'genes likelihood' file
		//---
		String entryId = tr.getChromosomeName() + ":" + start + "-" + end + "[" + aaExpected + "]";

		// Extract codon
		String cdsSeq = tr.cds();
		String codonStr = cdsSeq.substring(cdsBase, cdsBase + 3);
		String aa = tr.codonTable().aa(codonStr);

		if (aa.equals("" + aaExpected)) {
			if (debug) Gpr.debug("OK: " + entryId + " : " + aa);
		} else {
			if (debug) Gpr.debug("Entry ID     : " + entryId //
					+ "\ntr ID        : " + trId + ", chr: " + tr.getChromosomeName() + ", start: " + start + ", end: " + end + ", idx: " + aaIdx + ", fr: " + fr//
					+ "\nTranscript : " + tr //
					+ "\nExon       : " + ex //
					+ "\nStart pos: " + startPos //
					+ "\nCodon    : " + codonStr + ", aa (real): " + aa + ", aa (exp): " + aaExpected //
					);
			return false;
		}

		//---
		// Set marker coordinates
		//---
		parent = tr.getChromosome();
		if (tr.isStrandPlus()) {
			start = cds2pos[cdsBase];
			end = cds2pos[cdsBase + 2];
		} else {
			start = cds2pos[cdsBase + 2];
			end = cds2pos[cdsBase];
		}

		return true;
	}

	public void setMsa(String msaId, int aaIDx) {
		this.msaId = msaId;
		aaIdx = aaIDx;
	}

}
