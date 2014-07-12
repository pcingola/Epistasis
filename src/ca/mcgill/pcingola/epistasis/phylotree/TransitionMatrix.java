package ca.mcgill.pcingola.epistasis.phylotree;

import jeigen.DenseMatrix;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Calculate transition matrix.
 *
 * Based on two different libraries:
 * 		- Apache math commons
 * 		- JEigen: https://github.com/hughperkins/jeigen
 *
 * Note: We use JEigen for matrix exponential and log functions because Apache commons
 *       simply cannot do it (it cannot even decompose matrices with complex
 *       eigenvalues).
 *       I compared JEigen's results to R and they match within 10^-15 (at least
 *       for the simple 20x20 examples I used).
 *
 * @author pcingola
 */
public class TransitionMatrix extends Array2DRowRealMatrix {

	private static final long serialVersionUID = 1L;
	public static final double ACCEPTED_ERROR = 1e-4;
	public static final double EPSILON = 1e-6;

	protected EigenDecomposition eigen;
	protected boolean checkNegativeLambda;
	protected String colNames[];
	protected String rowNames[];

	/**
	 * Load from file
	 */
	public static double[][] loadValues(String fileName) {
		String file = Gpr.readFile(fileName);
		String lines[] = file.split("\n");

		int numRows = lines.length;
		int numCols = lines[0].split("\t").length;
		double d[][] = new double[numRows][numCols];

		int row = 0;
		for (String line : lines) {
			String recs[] = line.split("\t");
			for (int col = 0; col < recs.length; col++)
				d[row][col] = Gpr.parseDoubleSafe(recs[col]);
			row++;
		}

		return d;
	}

	public TransitionMatrix(double matrix[][]) {
		super(matrix);
	}

	public TransitionMatrix(int matrix[][]) {
		super(matrix.length, matrix[0].length);
		for (int i = 0; i < matrix.length; i++)
			for (int j = 0; j < matrix.length; j++)
				setEntry(i, j, matrix[i][j]);
	}

	public TransitionMatrix(int N) {
		super(new double[N][N]);
	}

	public TransitionMatrix(int rows, int cols) {
		super(rows, cols);
	}

	public TransitionMatrix(RealMatrix m) {
		super(m.getData());
	}

	public TransitionMatrix(String fileName) {
		super(loadValues(fileName));
	}

	public TransitionMatrix add(final TransitionMatrix m) throws MatrixDimensionMismatchException {
		// Safety check.
		MatrixUtils.checkAdditionCompatible(this, m);

		final int rowCount = getRowDimension();
		final int columnCount = getColumnDimension();
		final double[][] data = getData();
		final double[][] mdata = m.getData();
		final double[][] outData = new double[rowCount][columnCount];
		for (int row = 0; row < rowCount; row++) {
			final double[] dataRow = data[row];
			final double[] mRow = mdata[row];
			final double[] outDataRow = outData[row];
			for (int col = 0; col < columnCount; col++) {
				outDataRow[col] = dataRow[col] + mRow[col];
			}
		}

		return new TransitionMatrix(outData);
	}

	/**
	 * Show matrix's eigenvalues
	 */
	public double checkEien(boolean verbose) {
		eigen();

		// Exponentiate the diagonal
		RealMatrix D = eigen.getD().copy();
		double maxLambda = Double.NEGATIVE_INFINITY;
		int dim = D.getColumnDimension();
		for (int i = 0; i < dim; i++) {
			double lambda = D.getEntry(i, i);
			maxLambda = Math.max(maxLambda, lambda);
			if (verbose) System.out.println("\tlambda_" + i + "\t" + lambda + "\t" + eigen.getEigenvector(i));
		}

		if (verbose) System.out.println("\tlambda_max:\t" + maxLambda);
		if (maxLambda > 0) Gpr.debug("All Q's eigenvalues should be non-positive!");

		return maxLambda;
	}

	/**
	 * Perform eigen-decomposition and some sanity checks
	 */
	protected void eigen() {
		if (eigen == null) eigen = new EigenDecomposition(this);
	}

	/**
	 * Matrix exponentiation
	 */
	public RealMatrix exp(double time) {
		// Use Jeigen to calculate matrix log (using native methods)
		DenseMatrix m = new DenseMatrix(getData());
		DenseMatrix res = m.mexp();

		// Copy results from Jeigen
		int rows = getRowDimension();
		int cols = getColumnDimension();
		double d[][] = new double[rows][cols];
		for (int i = 0; i < getRowDimension(); i++)
			for (int j = 0; j < getColumnDimension(); j++)
				d[i][j] = res.get(i, j);

		return new TransitionMatrix(d);
	}

	public String[] getColNames() {
		return colNames;
	}

	public String[] getRowNames() {
		return rowNames;
	}

	/**
	 * Does this matrix have complex eigenvalues?
	 */
	public boolean hasComplexEigenvalues() {
		eigen();
		return eigen.hasComplexEigenvalues();
	}

	/**
	 * Does this matrix have any negative element off diagonal?
	 */
	public boolean hasNegativeOffDiagonalEntries() {
		int rows = getRowDimension();
		int cols = getColumnDimension();

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (i != j && getEntry(i, j) < 0) return true;
			}
		}
		return false;
	}

	/**
	 * Is this matrix diagonal?
	 */
	public boolean isDiagonal() {
		int rows = getRowDimension();
		int cols = getColumnDimension();

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				if (i != j && Math.abs(getEntry(i, j)) > EPSILON) return false;

		return true;
	}

	/**
	 * Is this matrix an identity matrix
	 */
	public boolean isIdentity() {
		int rows = getRowDimension();
		int cols = getColumnDimension();

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++) {
				if (i != j && Math.abs(getEntry(i, j)) > EPSILON) return false;
				if (i == j && Math.abs(getEntry(i, j) - 1.0) > EPSILON) return false;
			}

		return true;
	}

	/**
	 * Is this a probability matrix?
	 */
	public boolean isProbabilityMatrix() {
		int rows = getRowDimension();
		int cols = getColumnDimension();

		for (int i = 0; i < rows; i++) {
			double sum = 0;
			for (int j = 0; j < cols; j++) {
				double d = getEntry(i, j);
				if ((d < 0) || (d > 1)) return false;
				sum += d;
			}

			if (Math.abs(sum - 1) > EPSILON) return false;
		}

		return true;
	}

	/**
	 * Is this matrix symmetric?
	 */
	public boolean isSymmetric() {
		int rows = getRowDimension();
		int cols = getColumnDimension();

		for (int i = 0; i < rows; i++) {
			for (int j = i + 1; j < cols; j++) {
				double dij = getEntry(i, j);
				double dji = getEntry(j, i);
				double maxAbs = Math.max(Math.abs(dij), Math.abs(dji));

				if (maxAbs > 0) {
					double diff = Math.abs(dij - dji) / maxAbs;
					if (diff > EPSILON) return false;
				}
			}
		}

		return true;
	}

	/**
	 * Is this matrix zero?
	 */
	public boolean isZero() {
		int rows = getRowDimension();
		int cols = getColumnDimension();

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (getEntry(i, j) != 0) return false;
			}
		}

		return true;
	}

	/**
	 * Matrix log (natural log) times 1/time
	 */
	public RealMatrix log() {
		// Use Jeigen to calculate matrix log (using naitive methods)
		DenseMatrix m = new DenseMatrix(getData());
		DenseMatrix res = m.mlog();

		// Copy results from Jeigen
		int rows = getRowDimension();
		int cols = getColumnDimension();
		double d[][] = new double[rows][cols];
		for (int i = 0; i < getRowDimension(); i++)
			for (int j = 0; j < getColumnDimension(); j++)
				d[i][j] = res.get(i, j);

		return new TransitionMatrix(d);
	}

	public RealMatrix matrix(double time) {
		return this;
	}

	/**
	 * Save data to file
	 */
	public void save(String fileName) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < getRowDimension(); i++) {
			for (int j = 0; j < getColumnDimension(); j++) {
				if (j > 0) sb.append('\t');
				sb.append(getEntry(i, j));
			}
			sb.append('\n');
		}

		Gpr.toFile(fileName, sb.toString());
	}

	public void setColNames(String[] colNames) {
		this.colNames = colNames;
	}

	public void setNames(TransitionMatrix m) {
		colNames = m.colNames;
		rowNames = m.rowNames;
	}

	public void setRowNames(String[] rowNames) {
		this.rowNames = rowNames;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		// Column title
		if (colNames != null) {
			for (int i = 0; i < colNames.length; i++)
				sb.append((i > 0 ? "\t" : "") + colNames[i]);
			sb.append("\n");
		}

		// Data
		for (int i = 0; i < getRowDimension(); i++) {
			if (rowNames != null) sb.append(rowNames[i] + '\t');
			for (int j = 0; j < getColumnDimension(); j++) {
				if (j > 0) sb.append("\t");
				double val = getEntry(i, j);
				// sb.append(String.format("%1.6e", val));
				sb.append(val);
			}
			sb.append("\n");
		}

		return sb.toString();
	}

	public String toStringNice() {
		StringBuilder sb = new StringBuilder();

		// Column title
		if (colNames != null) {
			sb.append("        | ");
			for (int i = 0; i < colNames.length; i++)
				sb.append((i > 0 ? ", " : "") + String.format("%10s", colNames[i]));
			sb.append(" |\n");
		}

		// Data
		for (int i = 0; i < getRowDimension(); i++) {

			if (rowNames != null) sb.append(String.format("%8s| ", rowNames[i]));
			else sb.append("| ");

			for (int j = 0; j < getColumnDimension(); j++) {
				if (j > 0) sb.append(", ");
				double val = getEntry(i, j);
				double aval = Math.abs(val);
				if (aval < 1000000 && aval >= 100000.0) sb.append(String.format("% 6.2f", val));
				else if (aval < 100000 && aval >= 10000.0) sb.append(String.format("% 5.2f ", val));
				else if (aval < 10000 && aval >= 1000.0) sb.append(String.format("% 4.2f  ", val));
				else if (aval < 1000 && aval >= 100.0) sb.append(String.format("% 3.2f   ", val));
				else if (aval < 100 && aval >= 10.0) sb.append(String.format("% 2.2f    ", val));
				else if (aval < 10 && aval >= 1.0) sb.append(String.format("% 1.3f    ", val));
				else if (aval < 1.0 && aval >= 0.01) sb.append(String.format("% 1.3f    ", val));
				else if (aval < 1.0 && aval >= 0.001) sb.append(String.format("% 1.4f   ", val));
				else if (aval < 1.0 && aval >= 0.000001) sb.append(String.format("% 1.6f ", val));
				else if (aval < 1.0 && aval >= 0.0000001) sb.append(String.format("% 1.7f", val));
				else if (val == 0.0) sb.append(String.format(" 0        ", val));
				else sb.append(String.format("% 1.3e", val));
			}
			sb.append(" |\n");
		}

		return sb.toString();
	}

	/**
	 * Create a new (zero) matrix, same size as this one
	 */
	public TransitionMatrix zero() {
		int rows = getRowDimension();
		int cols = getColumnDimension();
		return new TransitionMatrix(rows, cols);
	}
}
