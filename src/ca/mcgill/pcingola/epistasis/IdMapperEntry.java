package ca.mcgill.pcingola.epistasis;

/**
 * An entry in a ID mapping file
 * 
 * @author pcingola
 */
public class IdMapperEntry {

	public String geneId, trId, geneName, refSeqId, pdbId;

	public IdMapperEntry(String line) {
		parseLine(line);
	}

	void parseLine(String line) {
		String fields[] = line.split("\t", -1);
		geneId = fields[0];
		trId = fields[1];
		geneName = fields[2];
		refSeqId = fields[3];
		pdbId = fields[4];
	}

	@Override
	public String toString() {
		return geneId + "\t" + trId + "\t" + geneName + "\t" + refSeqId + "\t" + pdbId;
	}

}
