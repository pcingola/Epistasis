package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.pcingola.epistasis.coordinates.GenomicCoordinates;

/**
 * Store genotype information
 *
 * @author pcingola
 */
public class Genotype extends GenomicCoordinates {

	private static final long serialVersionUID = 1L;

	protected int minorAlleleCount;
	protected byte gt[];

	public Genotype(Marker parent, int start, int end, String id, byte gt[]) {
		super(parent, start, end, id);
		this.gt = minorAllele(gt);
	}

	public Genotype(VcfEntry ve) {
		super(ve);
		gt = minorAllele(ve.getGenotypesScores());
		annotataions = ve.getInfo("EFF");
	}

	public byte[] getGt() {
		return gt;
	}

	public int getMinorAlleleCount() {
		return minorAlleleCount;
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
}
