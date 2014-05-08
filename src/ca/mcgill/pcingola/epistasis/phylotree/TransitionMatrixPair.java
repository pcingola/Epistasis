package ca.mcgill.pcingola.epistasis.phylotree;

import java.util.HashMap;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Transition matrix for amino acid pairs
 *
 * @author pcingola
 */
public class TransitionMatrixPair extends TransitionMatrixMarkov {

	public TransitionMatrixPair(double matrix[][]) {
		super(matrix);
	}

	public TransitionMatrixPair(RealMatrix m) {
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

}
