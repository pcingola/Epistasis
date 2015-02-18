package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.util.Gpr;
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

	public Genotype(Genome genome, String str) {
		super(null, 0, 0, str);
		parse(genome, str);
	}

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
		int countNonMIssing = 0;
		for (int i = 0; i < gt.length; i++)
			if (gt[i] >= 0) { // Skip missing alleles
				minorAlleleCount += gt[i]; // Allele count
				countNonMIssing++;
			}

		if (minorAlleleCount <= countNonMIssing) return gt; // OK, gt[] is minor allele

		// Convert to minor allele
		for (int i = 0; i < gt.length; i++)
			if (gt[i] >= 0) gt[i] = (byte) (2 - gt[i]);

		// Convert to minor allele
		minorAlleleCount = 2 * countNonMIssing - minorAlleleCount;

		return gt;
	}

	public int numberSamples() {
		return gt.length;
	}

	/**
	 * Parse a string formatted as '1:7724803_G/A'
	 */
	void parse(Genome genome, String str) {
		String f1[] = str.split(":");

		parent = genome.getOrCreateChromosome(f1[0]);

		String f2[] = f1[1].split("_");
		start = end = Gpr.parseIntSafe(f2[0]);
	}
}
