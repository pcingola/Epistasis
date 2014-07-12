package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.HashMap;

import org.apache.commons.math3.linear.RealMatrix;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Calculate transition matrix
 *
 * Use Markov model: P(t) = exp( t * Q )
 *
 * @author pcingola
 */
public class TransitionMatrixMarkov extends TransitionMatrix {

	private static final long serialVersionUID = 1L;

	HashMap<Double, RealMatrix> matrixByTime; // Cache matrices

	/***
	 * Load from file
	 */
	public static TransitionMatrixMarkov load(String fileName) {
		String file = Gpr.readFile(fileName);
		if (file == null || file.isEmpty()) throw new RuntimeException("Cannot read data from file '" + fileName + "'");

		String lines[] = file.split("\n");

		int numRows = lines.length - 1;
		int numCols = lines[0].split("\t").length;
		double d[][] = new double[numRows][numCols];
		String colNames[] = null;
		String rowNames[] = new String[numRows];

		int row = 0;
		for (String line : lines) {
			String recs[] = line.split("\t");

			if (row == 0) {
				// Title row
				colNames = recs;
			} else {
				for (int col = 0; col < recs.length; col++) {
					if (col == 0) rowNames[row - 1] = recs[col];
					else d[row - 1][col - 1] = Gpr.parseDoubleSafe(recs[col]);
				}
			}

			row++;
		}

		TransitionMatrixMarkov m = new TransitionMatrixMarkov(d);
		m.setColNames(colNames);
		m.setRowNames(rowNames);
		return m;
	}

	public TransitionMatrixMarkov(double matrix[][]) {
		super(matrix);
	}

	public TransitionMatrixMarkov(RealMatrix m) {
		super(m.getData());
	}

	@Override
	public RealMatrix matrix(double time) {
		// Check cache
		if (matrixByTime == null) matrixByTime = new HashMap<Double, RealMatrix>();
		RealMatrix m = matrixByTime.get(time);
		if (m != null) return m;

		// Perform matrix exponential
		m = exp(time);

		// Sanity check
		if (checkNegativeLambda) {
			double min = Double.POSITIVE_INFINITY;
			double max = Double.NEGATIVE_INFINITY;
			for (int i = 0; i < m.getRowDimension(); i++)
				for (int j = 0; j < m.getColumnDimension(); j++) {
					double e = m.getEntry(i, j);

					min = Math.min(min, e);
					max = Math.max(max, e);

					// Make sure we don't have numbers outside valid range
					if (e < 0) m.setEntry(i, j, 0.0);
					else if (e > 1.0) m.setEntry(i, j, 1.0);

				}

			if (min < -ACCEPTED_ERROR || max > (1.0 + ACCEPTED_ERROR)) throw new RuntimeException("All entries should be in [0, 1]\n\tmin : " + min + "\n\tmax : " + max);
		}

		// Add to cache
		matrixByTime.put(time, m);
		return m;
	}

	public void setCheck(boolean check) {
		checkNegativeLambda = check;
	}
}
