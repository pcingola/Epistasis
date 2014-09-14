package meshi.symmetryComplex.transformations;

/**
 * This class represents a transformation on points in space.
 * A 'legal' transformation can only be composed from translation and/or
 * rotation transformation, even though scale transformations exist.
 */
public class Transformation extends Matrix1D {

	/* Default transformation is no transformation using the unit matrix. */
	public Transformation() {
		super (4,4);
		
		for (int i=0; i<=3; i++)
			this.setMember(i, i, 1.0);
	}
	
	public Transformation(Matrix other) {
		super (other);

		if (! isLegalTransformation(other) )
			throw new RuntimeException
			("Illegal transformation:\n"+other);
	}

	/**
	 * Transformation constructor which can handle double[4][4]
	 * and double[3][3]. Everything else causes a RuntimeException,
	 * by MatrixAsDouble2D if the argument isn't a matrix at all
	 * else by this constructor.
	 * @param arrayMatrix the <code>double[][]</code> to be used for
	 * constructing this Transformation.
	 */
	public Transformation(double[][] arrayMatrix) {
		super(arrayMatrix);
		
//		if (! isLegalTransformation(this) )
//			throw new RuntimeException("Illegal matrix:\n"+this);

		if (rows==3 && (columns==3 | columns==4))
			expandToSquareFour(arrayMatrix);

/*
		if (rows==3 &&columns==4) {
			threeByFourToSquareFour(arrayMatrix);
		}
*/

		if (rows==4 && columns==4 && isLegalTransformation(this))
			return;

		throw new RuntimeException("parameter must be "+
				"a 4x4 or 3x3 or 3x4 legal transformation.\n");
	}

	/**
	 * Axis-angle constructor. The axis is given as xyz coordinates.
	 * Assumes that the magnitude of the axis is 1 and that theta is given 
	 * in radians. Code adapted from 
	 * <A href="http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/index.htm">here</A>.
	 */
	public Transformation(double x, double y, double z, double theta) {
		this();
		double c = Math.cos(theta);
		double s = Math.sin(theta);
		double t = 1.0 - c;
		
		setMember(0, 0, c + x*x*t);
		setMember(1, 1, c + y*y*t);
		setMember(2, 2, c + z*z*t);
		
		double tmp1 = x*y*t;
		double tmp2 = z*s;
		setMember(1, 0, tmp1 + tmp2);
		setMember(0, 1, tmp1 - tmp2);
		tmp1 = x*z*t;
		tmp2 = y*s;
		setMember(2, 0, tmp1 - tmp2);
		setMember(0, 2, tmp1 + tmp2);
		tmp1 = y*z*t;
		tmp2 = x*s;
		setMember(2, 1, tmp1 + tmp2);
		setMember(1, 2, tmp1 - tmp2);
	}

	private void expandToSquareFour(double[][] arrayMatrix) {
		boolean fourthRow = (arrayMatrix[0].length == 4);
		initFields(4, 4);
//		rows=4;
//		columns=4;
//		matrix = new double[4][4];
		for (int i=0; i<=2; i++) {
			for (int j=0; j<(fourthRow ? 4 : 3); j++)
				setMember(i, j, arrayMatrix[i][j]);
//			System.arraycopy(arrayMatrix[i], 0,
//					matrix[i], 0, fourthRow ? 4 : 3);
//			matrix[i][3] = fourthRow ? matrix[i][3] : 0.0;
		}
		setMember(3, 3, 1.0);
	}

/*
	private void threeByFourToSquareFour(double[][] arrayMatrix) {
		rows=4;
		columns=4;
		matrix = new double[4][4];
		for (int i=0; i<=2; i++) {
			System.arraycopy(arrayMatrix[i], 0, matrix[i], 0, 4);
			matrix[3][i] = 0.0;
		}
		matrix[3][3] = 1.0;
	}

	private void squareThreeToSquareFour(double[][] arrayMatrix) {
		rows=4;
		columns=4;
		matrix = new double[4][4];
		for (int i=0; i<=2; i++) {
			System.arraycopy(arrayMatrix[i], 0, matrix[i], 0, 3);
			matrix[i][3] = 0.0;
			matrix[3][i] = 0.0;
		}
		matrix[3][3] = 1.0;
	}
*/

	/**
	 * A legal transformation is composed of translation and/or
	 * rotation only.<BR>
	 * Such matrices have the following attribute:
	 * The 3x3 submatrix [0..2]x[0..2] is orthogonal.<BR>
	 * Orthogonality is tested in this way: A given matrix is
	 * orthogonal iff when multiplying it by its transposed matrix
	 * the result is the unit matrix.
	 */
	public static boolean isLegalTransformation(Matrix other) {
		return (other != null &&
				other.getRows()==4 &&
				other.getColumns()==4 &&
				other.getSubmatrix(0, 2, 0, 2).isOrthogonal());
	}
	
	public Transformation compose(Transformation other) {
		return new Transformation(other.multiply (this));
	}

	public double[][] transform(double[][] xyz) {
		if ( xyz==null || xyz.length!=4)
			throw new RuntimeException
			("Illegal array; must have exactly 4 members.");
		
		Matrix matrix = double2DToMatrix(xyz);
		matrix = this.multiply(matrix);
		
		return matrix.matrixToDouble2D();
	}
	
	public double[] transform(double[] xyz) {
		if (xyz==null || xyz.length!=3)
			throw new RuntimeException
			("Illegal array; must have exactly 3 members.");
		double[][] newXYZ = {{xyz[0]}, {xyz[1]}, {xyz[2]}, {1.0}};
		newXYZ = transform(newXYZ);
		return new double[] {newXYZ[0][0], newXYZ[1][0], newXYZ[2][0]};
	}


}