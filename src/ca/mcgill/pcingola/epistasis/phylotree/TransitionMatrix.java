package ca.mcgill.pcingola.epistasis.phylotree;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Calculate transition matrix
 * This trivial implementation simply returns exactly the same matrix for any 'time'
 *
 * @author pcingola
 */
public class TransitionMatrix extends Array2DRowRealMatrix {

	/**
	 * Load from file
	 */
	public static double[][] load(String fileName) {
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

	private static final long serialVersionUID = 1L;

	public static final double ACCEPTED_ERROR = 1e-4;
	protected EigenDecomposition eigen;
	protected boolean checkNegativeLambda;
	protected String colNames[];

	protected String rowNames[];

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
		super(load(fileName));
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
	public boolean checkEien(boolean verbose) {
		// Did we already perform eigendecomposition?
		if (eigen == null) eigen = new EigenDecomposition(this);

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
		if (maxLambda > 0) {
			Gpr.debug("All Q's eigenvalues should be non-positive!");
			return false;
		}

		return true;
	}

	/**
	 * Matrix exponentiation
	 */
	public RealMatrix exp(double time) {
		// Did we already perform Eigen-decomposition?
		if (eigen == null) eigen = new EigenDecomposition(this);

		// Exponentiate the diagonal
		RealMatrix D = eigen.getD().copy();
		int dim = D.getColumnDimension();
		RealMatrix expD = new DiagonalMatrix(dim);
		double maxLambda = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < dim; i++) {
			double lambda = D.getEntry(i, i);
			maxLambda = Math.max(maxLambda, lambda);
			expD.setEntry(i, i, Math.exp(lambda * time));
		}

		if (checkNegativeLambda && maxLambda > ACCEPTED_ERROR) throw new RuntimeException("All eigenvalues should be negative: max(lambda) = " + maxLambda);

		// Perform matrix exponential
		return eigen.getV().multiply(expD).multiply(eigen.getVT());
	}

	public String[] getColNames() {
		return colNames;
	}

	public String[] getRowNames() {
		return rowNames;
	}

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
		// Did we already perform eigendecomposition?
		if (eigen == null) eigen = new EigenDecomposition(this);

		// Exponentiate the diagonal
		RealMatrix D = eigen.getD();
		int dim = D.getColumnDimension();

		RealMatrix logD = new DiagonalMatrix(dim);
		double min = Double.MAX_VALUE;
		for (int i = 0; i < dim; i++) {
			double lambda = D.getEntry(i, i);
			if (lambda > 0) logD.setEntry(i, i, Math.log(lambda));
			min = Math.min(min, lambda);
			//			else {
			//				Gpr.debug("Negative eigenvalue when calculating log: lambda = " + lambda);
			//				return new TransitionMatrix(dim);
			//			}
		}

		// Strategy 1: Replace negative entries by 'lambdaMin'
		double lambdaMin = min / 2;
		for (int i = 0; i < dim; i++) {
			double lambda = D.getEntry(i, i);
			if (lambda < 0) logD.setEntry(i, i, Math.log(lambdaMin));
		}

		// Perform matrix exponential
		return eigen.getV().multiply(logD).multiply(eigen.getVT());
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

	public String toStringice() {
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

}
