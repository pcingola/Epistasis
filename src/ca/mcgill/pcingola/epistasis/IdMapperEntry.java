package ca.mcgill.pcingola.epistasis;

/**
 * An entry in a ID mapping file
 * 
 * @author pcingola
 */
public class IdMapperEntry {

	public String geneId, trId, protId, geneName, ucscId, pfamId, pdbId;

	public IdMapperEntry(String line) {
		parseLine(line);
	}

	void parseLine(String line) {
		String fields[] = line.split("\t", -1);
		geneId = fields[0];
		trId = fields[1];
		protId = fields[2];
		geneName = fields[3];
		ucscId = fields[4];
		pfamId = fields[5];
		pdbId = fields[6];
	}

	@Override
	public String toString() {
		return geneId + "\t" + trId + "\t" + protId + "\t" + geneName + "\t" + ucscId + "\t" + pfamId + "\t" + pdbId;
	}

}
