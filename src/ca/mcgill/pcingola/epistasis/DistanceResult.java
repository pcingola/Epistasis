package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.bio.structure.AminoAcid;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class DistanceResult {

	// Pdb information
	public String pdbId;
	public String pdbChainId;
	public int aaPos1, aaPos2;
	public char aa1, aa2;
	public double distance;

	// Genomic information
	public String transcriptId;
	public String chr1, chr2;
	public int chr1Num, chr2Num;
	public int pos1, pos2;
	public String annotations1, annotations2;

	// MSA information
	public String aaSeq1, aaSeq2;
	public String msa1, msa2;
	public int msaIdx1, msaIdx2;

	public DistanceResult() {
		pdbId = pdbChainId = transcriptId = aaSeq1 = aaSeq2 = chr1 = chr2 = annotations1 = annotations2 = "";
		aaPos1 = aaPos2 = pos1 = pos2 = chr1Num = chr2Num - 1;
		distance = -1;
	}

	public DistanceResult(AminoAcid aa1, AminoAcid aa2, double distance) {
		this();
		setAa1(aa1);
		setAa2(aa2);
		this.distance = distance;
	}

	public DistanceResult(String line) {
		this();

		// Parse line
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
			chr1Num = Chromosome.number(chr1);
		}

		if (fields.length > n) {
			String chrPos2 = fields[n++];
			String f[] = chrPos2.split(":");
			chr2 = f[0];
			pos2 = Gpr.parseIntSafe(f[1]);
			chr2Num = Chromosome.number(chr1);
		}

		if (fields.length > n) transcriptId = fields[n++];
		if (fields.length > n) msa1 = fields[n++];
		if (fields.length > n) msa2 = fields[n++];
		if (fields.length > n) msaIdx1 = Gpr.parseIntSafe(fields[n++]);
		if (fields.length > n) msaIdx2 = Gpr.parseIntSafe(fields[n++]);
		if (fields.length > n) aaSeq1 = fields[n++];
		if (fields.length > n) aaSeq2 = fields[n++];
		if (fields.length > n) annotations1 = fields[n++];
		if (fields.length > n) annotations2 = fields[n++];
	}

	/**
	 * Compare by genomic position
	 */
	public int compareByPos(DistanceResult d) {
		// Compare first position
		int comp = Chromosome.compare(chr1, d.chr1);
		if (comp != 0) return comp;

		comp = pos1 - d.pos1;
		if (comp != 0) return comp;

		// Compare second position
		comp = Chromosome.number(chr2) - Chromosome.number(d.chr2);
		if (comp != 0) return comp;

		comp = Chromosome.compare(chr2, d.chr2);
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

	/**
	 * Return amino acid pair (sorted)
	 */
	public String getAaPair() {
		return aa1 <= aa2 ? aa1 + "-" + aa2 : aa2 + "-" + aa1;
	}

	/**
	 * Return amino acid pair (sorted) + all combinations of annotations
	 */
	public List<String> getAaPairAnnotations() {
		boolean reversed = (aa1 > aa2);
		String aaPair = reversed ? aa2 + "-" + aa1 : aa1 + "-" + aa2;

		List<String> anns = new ArrayList<>();
		Arrays.stream(annotations1.split(";")) //
		.forEach( //
				ann1 -> Arrays.stream(annotations2.split(";")) //
				.forEach( //
						ann2 -> anns.add(aaPair + "\t" //
								+ (reversed ? ann2 + "\t" + ann1 : ann1 + "\t" + ann2) //
								) //
						) //
				);

		return anns;
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
				+ "\t" + (!chr1.isEmpty() ? chr1 + ":" + pos1 : "") //
				+ "\t" + (!chr2.isEmpty() ? chr2 + ":" + pos2 : "") //
				+ "\t" + transcriptId //
				+ "\t" + msa1 //
				+ "\t" + msa2 //
				+ "\t" + msaIdx1 //
				+ "\t" + msaIdx2 //
				+ "\t" + aaSeq1 //
				+ "\t" + aaSeq2 //
				+ "\t" + annotations1 //
				+ "\t" + annotations2 //
				;
	}

	/**
	 * Show genomic positions only
	 */
	public String toStringPos() {
		return "" //
				+ (chr1 != null ? "\t" + chr1 + ":" + pos1 : "") //
				+ (chr2 != null ? "\t" + chr2 + ":" + pos2 : "") //
				;
	}
}
