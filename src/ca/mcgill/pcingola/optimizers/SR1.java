package ca.mcgill.pcingola.optimizers;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.optimizers.exceptions.LineSearchException;
import ca.mcgill.pcingola.optimizers.exceptions.OptimizerException;

/**
 *
 **/

public class SR1 extends Minimizer {

	// Constant parameters
	public final int MAX_NUM_VARIABLES = 3000;
	public static final double DEFAULT_ALLOWED_MAX_H_FACTOR = 100;
	public static final int DEFAULT_MAX_NUM_KICK_STARTS = 3; // Don't change this number unless necessary
	public static final int DEFAULT_NUM_STEP_STEEPEST_DECENT = 50;
	public static final double DEFAULT_STEP_SIZE_EXPANTION = 1.1; // 2.0;

	protected SteepestDecent steepestDecent;
	protected WolfeConditionLineSearch lineSearchWolfe;
	protected int n; // number of variables
	int np; // (n+1) * n/2 - size of H
	protected double[] Binvk; // Inverse Hessian approximation at iteration k
	protected double[] pk; // p_k : The search direction
	protected double[] xk; // x_k : The coordinates at iteration k
	protected double[] gradNegk; // The (-) gradient at iteration k
	protected double[] sk; // s_k : The coordinates difference before the inverse Hessian update
	protected double[] yk; // y_k : The (-) gradients difference before the inverse Hessian update
	protected double[] ak; // A_k = Hk * yk
	protected double[] coordinates; // The position and gradients of the system
	protected double[] bufferCoordinates;
	protected int iterationNum; // Iterations counter

	// BFGS parameters
	protected double allowedMaxH;

	// Wolf conditions line search parameters
	protected double c1;
	protected double c2;
	protected double extendAlphaFactorWolfSearch;
	protected static int maxNumEvaluationsWolfSearch;

	// Steepest descent module parameters
	int numStepsSteepestDecent;
	double initStepSteepestDecent;
	double stepSizeReductionSteepestDecent;
	double stepSizeExpansionSteepestDecent;

	public static double abs(double a) {
		if (a < 0) return -1 * a;
		return a;
	}

	public SR1(Energy energy) {
		this(energy, DEFAULT_ALLOWED_MAX_H_FACTOR * energy.getTheta().length //
		, DEFAULT_MAX_NUM_KICK_STARTS //
				, WolfeConditionLineSearch.DEFAULT_C1 //
				, WolfeConditionLineSearch.DEFAULT_C2//
				, WolfeConditionLineSearch.DEFAULT_EXTENDED_ALPHA_FACTOR //
				, WolfeConditionLineSearch.DEFAULT_MAX_NUM_EVALUATIONS //
				, DEFAULT_NUM_STEP_STEEPEST_DECENT //
				, SimpleStepLength.DEFAULT_INITIAL_STEP_LENGTH //
				, SimpleStepLength.DEFAULT_STEP_SIZE_REDUCTION //
				, DEFAULT_STEP_SIZE_EXPANTION //
		);
	}

	public SR1(Energy energy //
			, double allowedMaxH, int maxNumKickStarts // General minimization parameters
			, double c1, double c2, double extendAlphaFactorWolfSearch, int maxNumEvaluationsWolfSearch // Parameters specific to the Wolf conditions line search
			, int numStepsSteepestDecent, double initStepSteepestDecent, double stepSizeReductionSteepestDecent, double stepSizeExpansionSteepestDecent // Steepest Decent parameters
	) {
		super(energy);
		setParameters(allowedMaxH, maxNumKickStarts, c1, c2, extendAlphaFactorWolfSearch, maxNumEvaluationsWolfSearch, numStepsSteepestDecent, initStepSteepestDecent, stepSizeReductionSteepestDecent, stepSizeExpansionSteepestDecent);
	}

	@Override
	protected void init() throws OptimizerException {
		steepestDecent = new SteepestDecent(energy(), numStepsSteepestDecent, initStepSteepestDecent, stepSizeReductionSteepestDecent, stepSizeExpansionSteepestDecent);
		lineSearchWolfe = new WolfeConditionLineSearch(energy(), c1, c2, extendAlphaFactorWolfSearch, maxNumEvaluationsWolfSearch);
		coordinates = energy().getTheta();
		n = coordinates.length;
		np = (n + 1) * n / 2;
		bufferCoordinates = new double[n];
		pk = new double[n];
		xk = new double[n];
		yk = new double[n];
		gradNegk = new double[n];
		sk = new double[n];
		ak = new double[n];
		iterationNum = 0;

		// Debug other algorithms?
		//		steepestDecent.setDebug(debug);
		//		lineSearchWolfe.setDebug(debug);

		// Starting the BFGS minimization by a few steepest descent steps, followed
		// by inverse Hessian initialization (to identity matrix), gradients and
		// position (X) initialization
		kickStart();
	}

	/**
	 * Initializing the inverse Hessian to the identity matrix
	 */
	protected void initHessian() {
		int i, j, k = 0;
		Binvk = new double[np];
		for (i = 0; i < n; i++) {
			Binvk[k] = 1;
			k++;
			for (j = 0; j < (n - i - 1); j++) {
				Binvk[k] = 0;
				k++;
			}
		}
	}

	/**
	 * Starting the BFGS minimization by a few steepest descent steps, followed by inverse Hessian initialization
	 */
	@Override
	protected void kickStart() throws OptimizerException {
		if (debug) Gpr.debug("A kick start has occurred in iteration:" + iterationNum + "\n");

		// Run steepest descent
		steepestDecent.run();

		// Initialize Hessian
		iterationNum += numStepsSteepestDecent;
		initHessian();

		// Initialize Wolfe search
		lineSearchWolfe.reset(steepestDecent.lastStepLength());

		// Update energy
		energy().evaluate();
		energy.copyTheta(xk);
		energy.getGradient(-1, gradNegk);
	}

	/**
	 * BFGS algorithm
	 */
	@Override
	protected boolean minimizationStep() throws OptimizerException {
		// if (debug) Gpr.debug(this);

		double curv = 0; // The curvature index
		double ykBinvYkCurv; // Yk*Hk*Yk
		double coefSkSkT; // A temporary result
		double maxBinvk; // The maximal entry in H (in term of magnitude)
		int i, j, k; // auxilary counters
		double tempAbs;

		//---
		// Step 1: Calculate pk = Binvk * (- grad[ f(xk) ] )
		//---
		energy.getGradient(-1, gradNegk);

		for (i = 0; i < n; i++) {
			pk[i] = 0;
			k = i;

			for (j = 0; j < i; j++) {
				pk[i] += Binvk[k] * gradNegk[j];
				k += (n - 1 - j);
			}

			for (j = i; j < n; j++) {
				pk[i] += Binvk[k] * gradNegk[j];
				k++;
			}
		}

		//---
		// Step 2: Perform line search
		//---
		try {
			energy.copyTheta(bufferCoordinates);
			lineSearchWolfe.findStepLength();
		} catch (LineSearchException lsEx) {

			if (debug) {
				// Return the energy coordinates to those before the line search
				System.err.println("Line seach failed");
				System.err.println("exception code =  " + lsEx.code);
				System.err.println("exception message = " + lsEx.getMessage());

				lsEx.printStackTrace();
			}

			energy.setTheta(bufferCoordinates);
			energy().evaluate();
			return false;
		}

		//---
		// Step 3: Calculate
		//          i) grad[ f( x_(k+1) ) ],
		//          ii) s_k = alpha_k * p_k
		//          iii) y_k = grad[ f(x_(k+1)) ] - grad[ f(x_k) ]
		//          iv) curvature = 1 / (s_k^T * y_k)
		// Check for pathological curvature
		//---
		curv = 0;
		double xk1[] = energy.getTheta(); // x_(k+1)
		double gradk1[] = energy.getGradient(); // grad[ f( x_(k+1) ) ]
		for (i = 0; i < n; i++) {
			// Note: We have to calculate yk = grad[ f(x_(k+1)) ] - grad[ f(x_k) ]
			// But gradNegk = - grad[ f(x_k) ], so we have to add
			yk[i] = gradk1[i] + gradNegk[i];

			// Calculate sk = alpha_k * p_k
			// Since x_(k+1) = x_k + alpha_k * p_k = x_k + s_k
			//       => s_k = x_(k+1) - x_k
			sk[i] = xk1[i] - xk[i];

			// Update gradient and x_k for next iteration
			gradNegk[i] = -gradk1[i];
			xk[i] = xk1[i];

			// Curvature = y_k^T * s_k = s_k^T * y_k
			curv += yk[i] * sk[i];
		}

		if (curv == 0) {
			if (debug) Gpr.debug("Minimization Error: The inverse Hessian is very badly scaled, and is unreliable, curv : " + curv);
			return false;
		} else curv = -1 / curv;

		//---
		// Step 4: Update inverse Hessian estimation
		//         B_{k+1)^-1 = B_k^(-1)
		//                      + ( s_k^T * y_k + y_k^T Binv_k * y_k ) * 1 / (( s_k^T * y_k )^2) * ( s_k * s_k^T )
		//                      - ( Binv_k * y_k * s_k^T + s_k * y_k^T * Binv_k ) * 1 / ( s_k^T * y_k )
		//         Note that
		//               i)   curv = 1 / ( s_k^T * y_k )		is a scalar
		//               ii)  y_k * s_k^T = (s_k * y_k^T)^T		is a rank one matrix
		//               iii) s_k * s_k^T						is another rank one matrix
		//				 iV)  a_k = curv * Binv_k * y_k			is a vector
		//---

		// a_k =  Binv_k * y_k / ( s_k^T * y_k ) = curv * Binv_k * y_k
		for (i = 0; i < n; i++) {
			ak[i] = 0;
			k = i;

			for (j = 0; j < i; j++) {
				ak[i] += Binvk[k] * yk[j];
				k += (n - 1 - j);
			}

			for (j = i; j < n; j++) {
				ak[i] += Binvk[k] * yk[j];
				k++;
			}

			ak[i] = curv * ak[i];
		}

		// Calculate ykBinvYkCurv = y_k' * Binv_k * y_k
		//                    = y_k' * ( curv * Binv_k * y_k ) / curv
		//                    = y_k' * a_k / curv
		ykBinvYkCurv = 0;
		maxBinvk = 0;
		for (i = 0; i < n; i++)
			ykBinvYkCurv += yk[i] * ak[i];

		// Coefficient multiplying s_k * s_k^T matrix:
		//
		// coefSkSkT = ( s_k^T * y_k + y_k' * Binv_k * y_k ) / ( s_k^T * y_k )^2
		//           = ( s_k^T * y_k ) / ( s_k^T * y_k )^2 + ( y_k' * Binv_k * y_k ) / ( s_k^T * y_k )^2
		//           = 1 / ( s_k^T * y_k ) + ( y_k' * Binv_k * y_k ) / ( s_k^T * y_k )^2
		//           = curv + ( ykBinvYkCurv ) / ( s_k^T * y_k )
		//           = curv + ykBinvYkCurv * curv
		coefSkSkT = curv * ykBinvYkCurv + curv;
		k = 0;

		// Update inverse Hessian estimation
		//         B_{k+1)^-1 = B_k^(-1)
		//                      + ( s_k^T * y_k + y_k^T Binv_k * y_k ) * 1 / (( s_k^T * y_k )^2)
		//                      - ( Binv_k * y_k * s_k^T + s_k * y_k^T * Binv_k ) * 1 / ( s_k^T * y_k )
		//
		//                    = B_k^(-1) + coefSkSkT * ( s_k * s_k^T ) - ( a_k * s_k^T + s_k * a_k^T )

		// Binv_(k+1) = Binv_k + (s_k * a_k' + a_k * s_k') + coefSkSkT * s_k * s_k'
		for (i = 0; i < n; i++) {
			for (j = i; j < n; j++) {
				Binvk[k] = Binvk[k] + coefSkSkT * sk[i] * sk[j] - (ak[j] * sk[i] + ak[i] * sk[j]);
				tempAbs = Binvk[k] * Binvk[k];
				if (tempAbs > maxBinvk) maxBinvk = tempAbs;
				k++;
			}
		}

		if (maxBinvk > allowedMaxH) {
			if (debug) Gpr.debug("Minimization Error: The inverse Hessian is very badly scaled, and is unreliable, maxBinvk : " + maxBinvk);
			return false;
		}

		return true;
	}

	/**
	 * Set initial parameters
	 */
	protected void setParameters(double allowedMaxH, int maxNumKickStarts // General
			, double c1, double c2, double extendAlphaFactorWolfSearch, int maxNumEvaluationsWolfSearch // Wolfe
			, int numStepsSteepestDecent, double initStepSteepestDecent, double stepSizeReductionSteepestDecent, double stepSizeExpansionSteepestDecent // Steepest descent
	) {
		this.allowedMaxH = allowedMaxH * allowedMaxH; // Doubling it so Math.abs is not needed in the comparison
		this.c1 = c1;
		this.c2 = c2;
		this.extendAlphaFactorWolfSearch = extendAlphaFactorWolfSearch;
		SR1.maxNumEvaluationsWolfSearch = maxNumEvaluationsWolfSearch;
		this.numStepsSteepestDecent = numStepsSteepestDecent;
		if (this.numStepsSteepestDecent < 1) this.numStepsSteepestDecent = 1;
		this.initStepSteepestDecent = initStepSteepestDecent;
		this.stepSizeReductionSteepestDecent = stepSizeReductionSteepestDecent;
		this.stepSizeExpansionSteepestDecent = stepSizeExpansionSteepestDecent;

		// Checking if the minimizing problem is not too large
		if (n > MAX_NUM_VARIABLES) throw new RuntimeException("\n\nThe number of variables to be minimized is greater than the maximum\n" + "this minimizer can handle. Use a minimizer for large-scale problems such as LBFGS\n");
	}
}
