package ca.mcgill.pcingola.epistasis.phylotree;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Simple matrix for calculating a trivial example with two times
 *
 * @author pcingola
 */
public class TransitionMatrix2Times extends TransitionMatrix {

	private static final long serialVersionUID = 1L;
	double timeThreshold;
	Array2DRowRealMatrix matrix1, matrix2;

	public TransitionMatrix2Times(double matrix1[][], double matrix2[][], double timeThreshold) {
		super(matrix1);
		this.matrix1 = new Array2DRowRealMatrix(matrix1);
		this.matrix2 = new Array2DRowRealMatrix(matrix2);
		this.timeThreshold = timeThreshold;
	}

	@Override
	public RealMatrix matrix(double time) {
		return time <= timeThreshold ? matrix1 : matrix2;
	}

	@Override
	public String toString() {
		return "Time <= " + timeThreshold + "\t" + matrix1 //
				+ "\nTime > " + timeThreshold + "\t" + matrix2;
	}
}
