package ca.mcgill.pcingola.regression;


/**
 * Weighted linear regression
 *
 * Source: http://www.codeproject.com/Articles/25335/An-Algorithm-for-Weighted-Linear-Regression
 *
 * Adapted by pcingola
 */
public class WeightedLinearRegression {

	double[][] V; // Least squares and var/covar matrix
	public double[] coefficients; // Coefficients
	public double[] stdErrCoeff; // Std Error of coefficients
	double RYSQ; // Multiple correlation coefficient
	double stdErr; // Standard deviation of errors
	double regressF; // Fisher F statistic for regression
	double[] yCalc; // Calculated values of Y
	double[] yResiduals; // Residual values of Y

	public double[] getCoefficients() {
		return coefficients;
	}

	public double getRegressF() {
		return regressF;
	}

	public double getRYSQ() {
		return RYSQ;
	}

	public double getStdErr() {
		return stdErr;
	}

	public double[] getStdErrCoeff() {
		return stdErrCoeff;
	}

	public double[] getyResiduals() {
		return yResiduals;
	}

	/**
	 * Perform regression
	 *  y[j]   = j-th observed data point
	 *  X[i,j] = j-th value of the i-th independent variable
	 *  w[j]   = j-th weight value
	 */
	public boolean regress(double[] y, double[][] X, double[] w) {

		int M = y.length; // M = Number of data points
		int N = X[0].length; // N = Number of linear terms
		int df = M - N; // Degrees of freedom
		yCalc = new double[M];
		yResiduals = new double[M];

		// If not enough data, don't attempt regression
		if (df < 1) { return false; }
		V = new double[N][N];
		coefficients = new double[N];
		stdErrCoeff = new double[N];
		double[] B = new double[N]; // Vector for LSQ

		// Clear the matrices to start out
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				V[i][j] = 0;

		// Form Least Squares Matrix
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				V[i][j] = 0;
				for (int k = 0; k < M; k++)
					V[i][j] = V[i][j] + w[k] * X[k][i] * X[k][j];
				//V[i][j] = V[i][j] + W[k] * X[i][k] * X[j][k];
			}
			B[i] = 0;
			for (int k = 0; k < M; k++)
				B[i] = B[i] + w[k] * X[k][i] * y[k];
			// B[i] = B[i] + W[k] * X[i][k] * Y[k];
		}

		// V now contains the raw least squares matrix
		if (!symmetricMatrixInvert(V)) return false;

		// V now contains the inverted least square matrix
		// Matrix multpily to get coefficients C = VB
		for (int i = 0; i < N; i++) {
			coefficients[i] = 0;
			for (int j = 0; j < N; j++)
				coefficients[i] = coefficients[i] + V[i][j] * B[j];
		}

		// Calculate statistics
		double TSS = 0;
		double RSS = 0;
		double YBAR = 0;
		double WSUM = 0;
		for (int k = 0; k < M; k++) {
			YBAR = YBAR + w[k] * y[k];
			WSUM = WSUM + w[k];
		}
		YBAR = YBAR / WSUM;
		for (int k = 0; k < M; k++) {
			yCalc[k] = 0;

			for (int i = 0; i < N; i++)
				yCalc[k] = yCalc[k] + coefficients[i] * X[k][i];
			//			Ycalc[k] = Ycalc[k] + C[i] * X[i][k];

			yResiduals[k] = yCalc[k] - y[k];
			TSS = TSS + w[k] * (y[k] - YBAR) * (y[k] - YBAR);
			RSS = RSS + w[k] * yResiduals[k] * yResiduals[k];
		}
		double SSQ = RSS / df;
		RYSQ = 1 - RSS / TSS;
		regressF = 9999999;
		if (RYSQ < 0.9999999) regressF = RYSQ / (1 - RYSQ) * df / (N - 1);
		stdErr = Math.sqrt(SSQ);

		// Calculate var-covar matrix and std error of coefficients
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++)
				V[i][j] = V[i][j] * SSQ;
			stdErrCoeff[i] = Math.sqrt(V[i][i]);
		}

		return true;
	}

	public boolean symmetricMatrixInvert(double[][] V) {
		int N = V.length;
		double[] t = new double[N];
		double[] Q = new double[N];
		double[] R = new double[N];
		double AB;
		int K, L, M;

		// Invert a symmetric matrix in V
		for (M = 0; M < N; M++)
			R[M] = 1;
		K = 0;

		for (M = 0; M < N; M++) {
			double big = 0;
			for (L = 0; L < N; L++) {
				AB = Math.abs(V[L][L]);
				if ((AB > big) && (R[L] != 0)) {
					big = AB;
					K = L;
				}
			}

			if (big == 0) return false;

			R[K] = 0;
			Q[K] = 1 / V[K][K];
			t[K] = 1;
			V[K][K] = 0;
			if (K != 0) {
				for (L = 0; L < K; L++) {
					t[L] = V[L][K];
					if (R[L] == 0) Q[L] = V[L][K] * Q[K];
					else Q[L] = -V[L][K] * Q[K];
					V[L][K] = 0;
				}
			}
			if ((K + 1) < N) {
				for (L = K + 1; L < N; L++) {
					if (R[L] != 0) t[L] = V[K][L];
					else t[L] = -V[K][L];
					Q[L] = -V[K][L] * Q[K];
					V[K][L] = 0;
				}
			}
			for (L = 0; L < N; L++)
				for (K = L; K < N; K++)
					V[L][K] = V[L][K] + t[L] * Q[K];
		}
		M = N;
		L = N - 1;
		for (K = 1; K < N; K++) {
			M = M - 1;
			L = L - 1;
			for (int J = 0; J <= L; J++)
				V[M][J] = V[J][M];
		}

		return true;
	}
}
