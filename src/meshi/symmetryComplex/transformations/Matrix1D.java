package meshi.symmetryComplex.transformations;

/**
 * An implementation of the <code>{@link meshi.util.mathTools.Jama.Matrix}</code> interface
 * as a one-dimesional array.
 * 
 * @author Oren
 *
 */
public class Matrix1D implements Matrix {

	private static final double EPSILON = 1e-4;

	protected int rows;
	protected int columns;
	protected double[] matrix;

	public Matrix1D() {
		initFields(1,1);
	}

	public Matrix1D(int rows, int columns) {
		if (rows >= 1 & columns >= 1)
			initFields(rows, columns);
		else
			throw new RuntimeException("matrix dimesnsions must be 1 or greater.");
	}

	public Matrix1D(Matrix other) {
		if (! isLegalMatrixObject(other) )
			throw new RuntimeException("Illegal Matrix Object.\n"+other);

		initFields(other.getRows(), other.getColumns());

		for (int row = 0; row < rows; row++)
			for (int column = 0; column < columns; column++)
				setMember(row, column, other.getMember(row, column));
	}

	public Matrix1D(double[][] matrix) {
//		if ( !isLegalMatrixObject(matrix) )
//			throw new RuntimeException
//			("Illegal matrix:\n" +matrix);
		
		initFields(matrix.length, matrix[0].length);
		
		for (int row = 0; row < rows; row++)
			for (int column = 0; column < columns; column++)
				setMember(row, column, matrix[row][column]);
		
//		if (! isOrthogonal() )
//			throw new RuntimeException("Illegal matrix.");
	}

	protected void initFields(int rows, int columns) {
		this.rows = rows;
		this.columns = columns;

		this.matrix = new double[rows * columns];
		for (int cell = 0; cell < this.matrix.length; cell++)
			this.matrix[cell] = 0.0;
	}

	public int getColumns() {return columns;}
	public int getRows() {return rows;}

	public double getMember(int row, int column) {
		if ( isLegalIndices(row, column) )
			return matrix[row * columns + column];
		else {
			throw new RuntimeException("Illegal indices for matrix.");
		}
	}

	public void setMember(int row, int column, double value) {
		if ( isLegalIndices(row, column) )
			matrix[row * columns + column] = value;
		else
			throw new RuntimeException("Illegal indices to matrix.");
	}

	public Matrix getSubmatrix(
			int fromRow, int toRow, int fromColumn, int toColumn) {

		if (! isLegalSubmatrixBoundaries(fromRow, toRow, fromColumn, toColumn) )
			throw new RuntimeException ("Illegal submatrix boundaries.");

		Matrix result = 
			new Matrix1D(toRow-fromRow+1, toColumn-fromColumn+1);

		for (int row = fromRow; row <= toRow; row++)
			for (int column = fromColumn; column <= toColumn; column++)
				result.setMember(row-fromRow, column-fromColumn, getMember(row, column));

		return result;
	}

	public boolean isLegalSubmatrixBoundaries(
			int fromRow, int toRow, int fromColumn, int toColumn) {

		return (0 <= fromRow && fromRow <= toRow && toRow < getRows() &&
				0 <= fromColumn && fromColumn <= toColumn && toColumn < getColumns());
	}

	public boolean isOrthogonal() {
		Matrix transposed = this.multiply(this.transpose());
		return getUnitMatrix(getRows()).equals(transposed);
	}

	public static Matrix getUnitMatrix(int size) {
		if ( size < 1 )
			throw new RuntimeException ("Illegal size: " +size);

		Matrix result = new Matrix1D(size, size);

		for (int i = 0; i < size; i++)
			result.setMember(i, i, 1.0);

		return result;

	}

	public boolean isZero(double verySmall) {
		return Math.abs(verySmall) <= EPSILON;
	}

	public Matrix double2DToMatrix(double[][] arrayMatrix) {
		int rows = arrayMatrix.length;
		int columns = arrayMatrix[0].length;

		Matrix result = new Matrix1D(rows, columns);

		for (int row = 0; row < rows; row++)
			for (int column = 0; column < columns; column++)
				result.setMember(row, column, arrayMatrix[row][column]);

		return result;
	}

	public double[][] matrixToDouble2D() {
		int rows = getRows ();
		int columns = getColumns();

		double[][] result = new double[rows][columns];

		for (int row = 0; row < rows; row++)
			for (int column = 0; column < columns; column++)
				result[row][column] = getMember(row, column);

		return result;
	}

	public Matrix multiply(Matrix other) {
		if ( this.getColumns() != other.getRows() )
			throw new RuntimeException
			("this.columns != other.rows");

		int k = this.getRows();
		int n = other.getColumns();
		int m = this.getColumns(); // == other.getRows()

		Matrix result = new Matrix1D(k, n);
		double sum;

		for (int i = 0; i < k; i++)
			for (int j = 0; j < n; j++) {
				sum = 0.0;

				for (int h = 0; h < m; h++)
					sum += ( this.getMember(i, h) * other.getMember(h, j) );

				result.setMember(i, j, sum);
			}

		return result;
	}

	public Matrix transpose() {
		int rows = getColumns();
		int columns = getRows();

		Matrix result = new Matrix1D(rows, columns);

		for (int i=0; i<rows; i++)
			for (int j=0; j<columns; j++)
				result.setMember(i, j, getMember(j, i));

		return result;
	}

	public boolean isLegalIndices(int row, int column) {
		return (row >= 0 &&
				row < rows &&
				column >= 0 &&
				column < columns );
	}
	public boolean equals(Object obj) {
		if (! (obj instanceof Matrix) )
			return false;

		return this.equals((Matrix)obj);
	}

	public boolean equals(Matrix other) {
		if ( (this.getRows() != other.getRows()) || (this.getColumns() != other.getColumns()))
			return false;

		for (int row = 0; row < this.getRows(); row++)
			for (int column = 0; column < this.getColumns(); column++)
				if (!( isZero(this.getMember(row, column) - other.getMember(row, column)) ))
					return false;

		return true;
	}

	public String toString() {
		String str = "";

		for (int i=0; i<getRows(); i++) {
			str += getMember(i, 0);

			for (int j=1; j<getColumns(); j++)
				str += "\t" + getMember(i, j);

			str += "\n";
		}

		return str;
	}

	public boolean isLegalMatrixObject(Object matrixObject) {
		if ( matrixObject == null )
			return false;
		
		if ( !(matrixObject instanceof Matrix1D) )
			return false;

		Matrix1D matrix = (Matrix1D) matrixObject;
		return matrix.rows * matrix.columns == matrix.matrix.length;
	}	

}
