package meshi.symmetryComplex.transformations;


public interface Matrix {
	
	public int getRows();
	public int getColumns();
	
	public double getMember(int row, int column);
	
	public void setMember(int row, int column, double value);
	
	public boolean isLegalIndices(int row, int column);
	
	public boolean isLegalMatrixObject(Object matrixObject);
	
	public Matrix multiply(Matrix other);
	
	public Matrix transpose();
	
	public Matrix double2DToMatrix(double[][] arrayMatrix);
	
	public double[][] matrixToDouble2D();
	
	public Matrix getSubmatrix(
		int fromRow, int toRow, int fromColumn, int toColumn);
	
	public boolean isLegalSubmatrixBoundaries(
		int fromRow, int toRow, int fromColumn, int toColumn);
	
	public boolean isOrthogonal();
	
	public boolean equals(Object obj);
	
	public boolean isZero(double verySmall);
	
	public String toString();
}