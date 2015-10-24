package ca.mcgill.pcingola.epistasis.coEvolutionMetrics;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Similarity matrix for Amino Acids
 *
 * @author pcingola
 */
public class AaSimilarityMatrix extends Array2DRowRealMatrix {

	private static final long serialVersionUID = 1L;
	protected static final int N = GprSeq.AMINO_ACIDS.length;

	public AaSimilarityMatrix() {
		super(N, N);
	}

	public AaSimilarityMatrix(String fileName) {
		super(N, N);
		load(fileName);
	}

	public void load(String fileName) {
		String file = Gpr.readFile(fileName);
		String lines[] = file.split("\n");

		String colNames[] = lines[0].split("\t");

		for (int i = 1; i < lines.length; i++) {
			String recs[] = lines[i].split("\t");
			int rowNum = GprSeq.aa2Code(recs[0].charAt(0));

			for (int j = 1; j < recs.length; j++) {
				String colName = colNames[j];
				int colNum = GprSeq.aa2Code(colName.charAt(0));
				setEntry(rowNum, colNum, Gpr.parseDoubleSafe(recs[j]));
			}
		}
	}

	@Override
	public String toString() {
		return Gpr.toString(getData());
	}

}
