package ca.mcgill.pcingola.epistasis;

import java.util.function.Function;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An entry in a ID mapping file
 *
 * @author pcingola
 */
public class IdMapperEntry implements Cloneable {

	// Select ID function
	public static final Function<IdMapperEntry, String> IDME_TO_REFSEQ = ime -> ime.refSeqId;
	public static final Function<IdMapperEntry, String> IDME_TO_ENSEMBLID = ime -> ime.trId;

	public String geneId, trId, geneName, refSeqId, pdbId, pdbChainId;
	public int pdbAaLen, trAaLen;

	public IdMapperEntry(String line) {
		geneId = trId = geneName = refSeqId = pdbId = pdbChainId = "";
		pdbAaLen = trAaLen = 0;
		parseLine(line);
	}

	/**
	 * Clone the object and return a copy having a different chain ID
	 */
	public IdMapperEntry cloneAndSet(String chainId, int pdbAaLen, int trAaLen) {
		try {
			IdMapperEntry cloned = (IdMapperEntry) this.clone();
			cloned.pdbChainId = chainId;
			cloned.pdbAaLen = pdbAaLen;
			cloned.trAaLen = trAaLen;
			return cloned;
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Parse line from a file
	 */
	void parseLine(String line) {
		String fields[] = line.split("\t", -1);
		int n = 0;
		geneId = fields[n++];
		trId = fields[n++];
		geneName = fields[n++];
		refSeqId = fields[n++];
		pdbId = fields[n++];
		if (fields.length > n) pdbChainId = fields[n++];
		if (fields.length > n) pdbAaLen = Gpr.parseIntSafe(fields[n++]);
		if (fields.length > n) trAaLen = Gpr.parseIntSafe(fields[n++]);
	}

	@Override
	public String toString() {
		return geneId //
				+ "\t" + trId //
				+ "\t" + geneName //
				+ "\t" + refSeqId //
				+ "\t" + pdbId //
				+ "\t" + pdbChainId //
				+ "\t" + pdbAaLen //
				+ "\t" + trAaLen //
		;
	}

}
