package ca.mcgill.pcingola.epistasis;

import java.util.function.Function;

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

	public IdMapperEntry(String line) {
		parseLine(line);
	}

	/**
	 * Clone the object and return a copy having a different chain ID
	 * @param chainId
	 * @return
	 */
	public IdMapperEntry cloneAndSetChainId(String chainId) {
		try {
			IdMapperEntry cloned = (IdMapperEntry) this.clone();
			cloned.pdbChainId = chainId;
			return cloned;
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}

	void parseLine(String line) {
		String fields[] = line.split("\t", -1);
		geneId = fields[0];
		trId = fields[1];
		geneName = fields[2];
		refSeqId = fields[3];
		pdbId = fields[4];
		if (fields.length >= 6) pdbChainId = fields[5];
	}

	@Override
	public String toString() {
		return geneId + "\t" + trId + "\t" + geneName + "\t" + refSeqId + "\t" + pdbId + "\t" + pdbChainId;
	}

}
