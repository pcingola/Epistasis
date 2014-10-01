package meshi.optimizers;

import meshi.optimizers.exceptions.LineSearchException;
import meshi.optimizers.exceptions.OptimizerException;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 *This class implements a BFGS minimizer according to the scheme in: Numerical Optimization by J. Nocendal &
 *S. J. Wright, Springer 1999, pp 193-201.
 *
 * This class was written by Nir Kalisman as part of the MESHI package Kalisman et al. (2005) Bioinformatics 21:3931-3932
 *
 *The BFGS algorithm (general)
 *----------------------------
 *In Newton minimizers an approximation to the Hessian of the energy function at position Xk is calculated. Then finding the
 *inverse of that Hessian (Hk), and solving the equation Pk = -Hk*grad(Xk) gives a good search direction Pk. Later, a
 *line search procedure has to determine just how much to go in that direction (producing the scalar alpha_k). The new
 *position is given by: Xk+1 = Xk + alpha_k*Pk.In the BFGS method the inverse Hessian is not computed explicitly. Instead
 *it is updated each step by the values of Pk and the new gradient. The updating formula is (t - transpose):
 *Hk+1 = (I  - Rk*Sk*Ykt)Hk(I  - Rk*Yk*Skt) + Rk*Skt*Sk
 *where:
 *Sk = Xk+1 - Xk
 *Yk = grad(Xk+1)-grad(Xk)
 *Rk = 1/(Ykt*Sk)
 *
 *
 *The BFGS algorithm (specific implementation)
 *--------------------------------------------
 *1)To run this minimizer:
 *a) Instantiate this class with the desired minimization parameters.
 *b) Put the initial coordinates in the 'coordinates' variable at the 'energy' class.
 *c) Activate BFGS.run().
 *d) Check for thrown errors to see if the minimization succeeded.
 *e) The minimized position is in the 'coordinates' variable at the 'energy' class.
 *
 *2)This implementation creates a matrix (of doubles) size is 0.5*(n^2). Where n is the number of variables to minimize. If n is large
 *the memory load might be very great.
 *
 *3)The inverse Hessian matrix (H) which is symmetric is stored as a linear vector in the following way (to save space):
 *H(1,1:n) followed by H(2,2:n) followed by H(3,3:n) and so on until H(n,n).
 *
 *4)The BFGS algorithm is generally robust and efficient. With certain energy functions it might, however, become unstable.
 *The algorithm checks for instabilities during the run (the thresholds to some of the instabilities are given by the parameters
 *in the constructor). If such instability is discovered, the algorithm try to recover by changing to steepest descent
 *algorithm for a few steps (kick-start). If this option is activated too much it is a sign of a deeper problem, and an
 *informative error message is thrown.
 *
 *5)Good default values are given after the parameter name.
 *
 *6)The initial guess to the Hessian is the unity matrix. Better guesses are possible (see the reference).
 *
 *General minimization parameters
 *-------------------------------
 *- energy - pointer to an TotalEnergy object, where the energy function is.
 *- tolerance - 1e-6 - Minimization stops when the magnitude of the maximal gradient component drops below tolerance.
 *- maxSteps - 1000 - The maximal number of iteration steps allowed
 *
 *
 *Parameters Specific to the BFGS algorithm
 *-----------------------------------------
 *- allowedMaxH - (10-100)*n - In some energy function scenarios the inverse Hessian approximation might become unstable and
 *                       unrealiable, by having huge numbers in the H entries. This parameter sets a upper limit on the
 *                       matrix H entries. Higher values would lead to a new kick-start. This value should be somewhere
 *                       in the range (10-100)*n (lower is more conservative).
 *- maxNumKickStarts - 3 - If the minimzer become unstable for some reason, it could be restarted from the current position.
 *                         This parameter determined how many times this could happen before the minimization is aborted.
 *
 *
 *Parameters specific to the Wolf conditions line search
 *------------------------------------------------------
 *The BFGS algorithm requires a step length finder who finds a step length that also satisfies the Wolf conditions. See the
 *help of this specific line search for explaination of what these conditions are. The parameters of this line search are:
 *
 *- c1,c2 - 1e-4,0.9 - The two parameters of the Wolf conditions. Must satisfy: 0<c1<c2<1
 *- maxNumEvaluations - 10 - The maximal number of step length trails allowed. This gives an upper limit on the total number of
 *                        evaluations both in the bracketing and Zoom. If line search fails because this number was
 *                        exceeded try to enlarge 'extendAlphaFactor' or change the initial alpha guess. This error might
 *                        also be caused if the tolerance of the minimizer is set to be extremely small.
 *- extendAlphaFactor - 3 - After a certain step length fails, the next step length is taken as this multiple of the failing
 *                        last alpha.
 *
 *
 *Steepest Decent module
 *-----------------------
 *In two cases steepest descent minimization is done instead of BFGS.
 *1) All runs start with a certain number of steepest descent steps, because difficult scenarios for BFGS minimization
 *might occur at the start due to atom clashes.
 *2) If the normal operation of the minimizer is disturbed for some reason  (failing to produce a descent direction,
 *failing to satisfies the wolf conditions, etc.) another set of steepest descent steps (with similar parameters to
 *case 1) is attempted. If the normal operation is disturbed too many times, the minimization is aborted because
 *this is indicative of a more severe fault, most likely in the energy function.
 *
 *The steepest descent parameters are as follow:
 *- numSteepestDecent - 50 - The number of steepest descent steps to be taken. If this number is smaller than 1, than at
 *                      least one steepest descent step is done.
 *- initialStepLength - 1 - parameter of the steepest descent line search. The first step length to be tried after the
 *                      calculation of the first gradient. This parameter should normally be 1 unless very large gradients
 *						(such as clashing of VDW atoms) are expected in the first steps. In that case it should be  set to
 *                      a much smaller value (1e-4 or less).
 *- stepSizeReduction - 0.5 - parameter of the line search. The step length is multiplied by this factor if no reduction
 *                      in energy is achieved.
 *- stepSizeExpansion - 2 - parameter of the line search. The first step length tried is the step length from previous
 *                      line search multiplied by this factor. (Note that non-positive values to this parameter cause
 *                      special options to be called (see the SimpleStepLength class help).
 *
 *
 *Constant Parameters
 *-------------------
 *MAX_NUM_VARIABLES - 1000 - The maximal number of variables of the energy functions that can be minimized using this
 *minimizer. Since the BFGS minimizer uses a matrix with size 0.5*n^2 (n - number of variables), it is obvious that
 *too many variables will lead to insufficient memory problems. Therefore the minimization is aborted with an error
 *if the number of variables exceeds this limit.
 *
 *
 *Reference: http://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm
 *
 **/

public class BFGS extends Minimizer {

	// Constant parameters
	public final int MAX_NUM_VARIABLES = 3000;
	public static final double DEFAULT_ALLOWED_MAX_H_FACTOR = 100;
	public static final int DEFAULT_MAX_NUM_KICK_STARTS = 3; // Don't change this number unless necessary
	public static final int DEFAULT_NUM_STEP_STEEPEST_DECENT = 50;
	public static final double DEFAULT_STEP_SIZE_EXPANTION = 2.0;

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
	protected double[] Ak; // A_k = Hk * yk
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

	public BFGS(Energy energy) {
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

	public BFGS(Energy energy //
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
		Ak = new double[n];
		iterationNum = 0;

		steepestDecent.setDebug(debug);
		lineSearchWolfe.setDebug(debug);

		// Starting the BFGS minimization by a few steepest descent steps, followed by inverse Hessian, gradients (G),
		// and position (X) initialization
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
		if (debug) Gpr.debug("Steepest descent, step: " + steepestDecent.lastStepLength() + ", " + Gpr.toString(energy.getTheta()));

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
		if (debug) Gpr.debug(this);

		double curv = 0; // The curvature index
		double YHY; // Yk*Hk*Yk
		double Coef; // A temporary result
		double MaxH; // The maximal entry in H (in term of magnitude)
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
			System.err.println("Minimization Error: The inverse Hessian is very badly scaled, and is unreliable\n");
			return false;
		} else curv = -1 / curv;

		//---
		// Step 4: Updating the inverse Hessian estimation
		//         B_{k+1)^-1 = B_k^(-1)
		//                      + ( s_k^T * y_k + y_k^T Binv_k * y_k ) * 1 / (( s_k^T * y_k )^2)
		//                      - ( Binv_k * y_k * s_k^T + s_k * y_k^T * Binv_k ) * 1 / ( s_k^T * y_k )
		//         Note that
		//               i)   curv = 1 / ( s_k^T * y_k )		is a scalar
		//               ii)  y_k * s_k^T = (s_k * y_k^T)^T		is a rank one matrix
		//               iii) s_k * s_k^T						is another rank one matrix
		//				 iV)  A_k = curv * Binv_k * y_k			is a vector
		//---

		// A_k = curv * Binv_k * y_k = Binv_k * y_k / ( s_k^T * y_k )
		for (i = 0; i < n; i++) {
			Ak[i] = 0;
			k = i;

			for (j = 0; j < i; j++) {
				Ak[i] += Binvk[k] * yk[j];
				k += (n - 1 - j);
			}

			for (j = i; j < n; j++) {
				Ak[i] += Binvk[k] * yk[j];
				k++;
			}

			Ak[i] = curv * Ak[i];
		}

		// Hk+1 = Hk + (Sk*Ak' + Ak*Sk') + (Curv*(Yk'*Ak) + Curv)*Sk*Sk'
		YHY = 0;
		MaxH = 0;
		for (i = 0; i < n; i++)
			YHY += yk[i] * Ak[i];
		Coef = curv * YHY + curv;
		k = 0;
		for (i = 0; i < n; i++) {
			for (j = i; j < n; j++) {
				Binvk[k] = Binvk[k] + Ak[j] * sk[i] + Ak[i] * sk[j] + Coef * sk[i] * sk[j];
				tempAbs = Binvk[k] * Binvk[k];
				if (tempAbs > MaxH) MaxH = tempAbs;
				k++;
			}
		}
		if (MaxH > allowedMaxH) {
			System.out.println("Minimization Error: The inverse Hessian is very badly scaled, and is unreliable\n");
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
		BFGS.maxNumEvaluationsWolfSearch = maxNumEvaluationsWolfSearch;
		this.numStepsSteepestDecent = numStepsSteepestDecent;
		if (this.numStepsSteepestDecent < 1) this.numStepsSteepestDecent = 1;
		this.initStepSteepestDecent = initStepSteepestDecent;
		this.stepSizeReductionSteepestDecent = stepSizeReductionSteepestDecent;
		this.stepSizeExpansionSteepestDecent = stepSizeExpansionSteepestDecent;

		// Checking if the minimizing problem is not too large
		if (n > MAX_NUM_VARIABLES) throw new RuntimeException("\n\nThe number of variables to be minimized is greater than the maximum\n" + "this minimizer can handle. Use a minimizer for large-scale problems such as LBFGS\n");
	}
}
