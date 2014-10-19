package ca.mcgill.pcingola.optimizers;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.pcingola.optimizers.exceptions.LineSearchException;

/**
 *This class implements a line search that satisfies the Wolf conditions according to the scheme in: Numerical Optimization
 *by J. Nocendal & S. J. Wright, Springer 1999, pp 55-62.
 *
 *Being in position Xk, and having a descent search direction Pk, a step length (alpha) is sought, that leads to a new position Xk+1,
 *so that Xk+1 = Xk + alpha*Pk, and stisfies the 2 Wolf conditons:
 *1) E(Xk+1) < E(Xk) + c1*alpha*grad(Xk)*Pk   (sufficient decrease)
 *2) |grad(Xk+1)*Pk| < c2*|grad(Xk)*Pk|   (curvature condition)
 *The implementation follows the one in the above reference. In short, the variable alpha is iteratively enlarged, until it
 *brackets (with the previous iteration value) an interval that contain some points with the Wolf conditions. The function
 *Zoom is then called, that finds the exact alpha where the Wolf conditions hold. The initial choice of alpha (at the first
 *iteration) is the popular min{1 , -2*(E(Xk-1) - E(Xk))/(grad(Xk)*Pk)}. Other choices are possible to code if this one
 *leads to failure of the line search algorithm. This algorithm is safe guarded against infinite loop, by a condition on the
 *maximal number of energy evaluations. See bellow what to do if this condition is violated.
 *
 *How to use this class:
 *----------------------
 *a) Instantiate this class with the desired parameters.
 *b) In a position Xk: evaluate the energy gradients and coordinates to position Xk.
 *c) Run 'findStepLength(Vec[n][2])' where the first column in Vec is the position Xk, and the second
 *column is the direction Pk. This method returns the found step length. This method also changes the coordinates
 *in class 'energy' to the coordinates in the new (minimized) position: Xk+1 = Xk + (step_length)*Pk. Also the
 *gradients in the energy class are updated.
 *d) Check for thrown exceptions to make sure that the step length is correct.
 *
 *
 *
 *Parameters (Good default values are given after the parameter name)
 *----------
 *c1,c2 - 1e-4,0.9 - The two paramters of the Wolf conditions. Must satisfy: 0<c1<c2<1
 *maxNumEvaluations - 10 - The maximal number of step length trails allowed. This gives an upper limit on the total number of
 *                        evaluations both in the bracketing and Zoom. If line search fails because this number was
 *                        exceeded try to enlarge 'extendAlphaFactor' or change the initial alpha guess. This error might
 *                        also be caused if the tolerance of the minimizer is set to be extremely small.
 *extendAlphaFactor - 3 - After a certain step length fails, the next step length is taken as this multiple of the failing
 *                        last alpha.
 *
 *Constant Parameters
 *-------------------
 *interSafeGuardFactor - 20 - This factor is used in the calculation of the cubic interpolation safeguard. It's exact value
 *                            is not crucial to the minimizer performance. 20 is generally good.
 *
 *
 **/

public class WolfeConditionLineSearch extends LineSearch {

	public static final double DEFAULT_C1 = 1e-4;
	public static final double DEFAULT_C2 = 0.9;
	public static final double DEFAULT_EXTENDED_ALPHA_FACTOR = 3.0;
	public static final int DEFAULT_MAX_NUM_EVALUATIONS = 10;

	private double maxNumEvaluations;
	private double extendAlphaFactor;
	private double c1, c2;
	private final double interSafeGuardFactor = 20;
	private double interSafeGuard;
	private double energyOld, energyNew;
	private double[] x;
	private double[] xCopy;
	private int firstRun;
	private double alphaPrev, alpha, alphaFinal;
	private int n; // Number of variables
	private double e0, ePrev, e; // Energy values along Pk
	private double grad0, pkGradPrev, pkGrad; // Gradients along Pk
	private int numAlphaEvaluations;
	boolean stop = false;
	private int i;

	public WolfeConditionLineSearch(Energy energy) {
		this(energy, DEFAULT_C1, DEFAULT_C2, DEFAULT_EXTENDED_ALPHA_FACTOR, DEFAULT_MAX_NUM_EVALUATIONS);
	}

	public WolfeConditionLineSearch(Energy energy, double c1, double c2, double extendAlphaFactor, int maxNumEvaluations) {
		super(energy);
		this.c1 = c1;
		this.c2 = c2;
		this.maxNumEvaluations = maxNumEvaluations;
		this.extendAlphaFactor = extendAlphaFactor;
		interSafeGuard = extendAlphaFactor / interSafeGuardFactor;

		x = energy.getTheta(); // Local copy of coordinates
		n = x.length;
		xCopy = new double[n];

		reset();
	}

	@Override
	public double findStepLength() throws LineSearchException {
		energy.copyTheta(xCopy);
		return findStepLength(xCopy, energy.getGradient(-1.0));
	}

	/**
	 * Find optimization step satisfiying wolfe conditions
	 * @param x0 : Initial coordinates
	 * @param pk : Line of descent (e.g. negative gradient)
	 */
	public double findStepLength(double x0[], double pk[]) throws LineSearchException {
		if (debug) Gpr.debug(this + "\n\tpk : " + Gpr.toString(pk));

		// Initializing the run
		energyNew = energy.getEnergy();

		// Calculating  grad_0 = - p_k^T * grad[ f(x_k) ]
		grad0 = 0;
		double grad[] = energy.getGradient();
		for (i = 0; i < n; i++)
			grad0 += pk[i] * grad[i];

		// Sanity check: Is this in the descent direction?
		if (grad0 >= 0) {
			if (debug) Gpr.debug("Error: Gradient direction issue?");
			throw new LineSearchException(LineSearchException.NOT_A_DESCENT_DIRECTION, "\n\nThe search direction is not a descent direction. This problem might be caused by incorrect diffrentiation of the energy function or by numerical instabilities of the minimizing techniques (such as not fullfilling the Wolf condtions in BFGS).\n");
		}

		numAlphaEvaluations = -1;
		stop = false;
		e0 = energyNew;
		ePrev = e0;
		alphaPrev = 0;

		// Choosing the initial alpha guess
		if (firstRun != 0) {
			firstRun = 0;
		} else {
			alpha = 2.02 * (energyNew - energyOld) / grad0;
			if (alpha > 1) alpha = 1;
		}

		energyOld = energyNew; // For the next time line search is called

		// Bracketing the Wolf area
		while ((numAlphaEvaluations <= maxNumEvaluations) && (!stop)) {
			numAlphaEvaluations++;

			// Calculating left hand side of first Wolfe condition: eI1 = f(x0 + alpha * pk)
			for (i = 0; i < n; i++)
				x[i] = x0[i] + alpha * pk[i];
			energy.needsUpdate(); // Force energy function update
			e = energy.evaluate();

			// Calculating left hand side of second Wolfe condition: gradI1 = pk * grad[ f(x0 + alpha * pk) ]
			pkGrad = 0;
			grad = energy.getGradient();
			for (i = 0; i < n; i++)
				pkGrad += pk[i] * grad[i];

			if ((e > (e0 + c1 * alpha * grad0)) || ((e >= ePrev) && (numAlphaEvaluations > 0))) {
				if (debug) Gpr.debug("Wolfe condition I not satisfied: " //
						+ "\n\tf(x_k + alpha * p_k) <= f(x_k) + c1 * alpha * p_k * grad[ f(x_k) ]: " + (e <= (e0 + c1 * alpha * grad0)) //
						+ "\n\tx_k                                        :" + Gpr.toString(x0) //
						+ "\n\tp_k                                        :" + Gpr.toString(pk) //
						+ "\n\tp_k * grad[ f(x_k) ]                       :" + grad0 //
						+ "\n\talpha                                      :" + alpha //
						+ "\n\tx_k + alpha * p_k                          :" + Gpr.toString(x) //
						+ "\n\tf( x_k )                                   :" + e0 //
						+ "\n\tf( x_k + alpha * p_k )                     :" + e //
						+ "\n\tf(x_k) + c1 * alpha * p_k * grad[ f(x_k) ] :" + (e0 + c1 * alpha * grad0) //
				);
				zoom(x0, pk, false); // Calling the regular zoom
			} else {
				if (Math.abs(pkGrad) <= (-c2 * grad0)) {
					alphaFinal = alpha;
					stop = true;
				} else if (pkGrad >= 0) zoom(x0, pk, true); // Inverse Zoom
			}

			alphaPrev = alpha;
			pkGradPrev = pkGrad;
			ePrev = e;
			alpha = alpha * extendAlphaFactor;
		}

		if (numAlphaEvaluations > maxNumEvaluations) {
			if (debug) Gpr.debug("Error: Too many evaluations");

			energy.setTheta(x0); // Returning the coordinates to the original state
			energy.evaluate();
			throw new LineSearchException(LineSearchException.WOLF_CONDITION_NOT_MET, "\n\nWolf conditions not met. The line search did not converge, and exceeded the maximal number of step extensions allowed.\n");
		}

		if (debug) Gpr.debug("alpha_final: " + alphaFinal + "\t" + this);

		return alphaFinal;
	}

	public void reset() {
		firstRun = 1;
		alpha = 1.0;
	}

	public void reset(double resAlpha) {
		firstRun = 1;
		alpha = resAlpha;
	}

	// The function Zoom finds a step length satisfying the Wolf conditions, given the bracketing of alphaI and alphaI1.
	// It was separated into a different function to make the code more readable.
	private void zoom(double[] x0, double[] pk, boolean inv) {
		double alphaHi, alphaLow, alphaNew = 0;
		double eHi, eLow, eNew = 0;
		double gradHi, gradLow, gradNew = 0;
		double a = 0, b = 0, ga = 0, gb = 0, ea, eb, interval, d1 = 0, d2 = 0; // Cubic interpolation variables

		eHi = e;
		gradHi = pkGrad;
		alphaHi = alpha;
		eLow = ePrev;
		gradLow = pkGradPrev;
		alphaLow = alphaPrev;

		if (inv) {
			eHi = ePrev;
			gradHi = pkGradPrev;
			alphaHi = alphaPrev;
			eLow = e;
			gradLow = pkGrad;
			alphaLow = alpha;
		}

		while ((numAlphaEvaluations <= maxNumEvaluations) && (!stop)) {
			numAlphaEvaluations++;

			// Cubic interpolation of the next step length
			a = alphaLow;
			b = alphaHi;
			ga = gradLow;
			gb = gradHi;
			ea = eLow;
			eb = eHi;

			if (a > b) { //switch
				b = alphaLow;
				a = alphaHi;
				gb = gradLow;
				ga = gradHi;
				eb = eLow;
				ea = eHi;
			}

			d1 = ga + gb - 3 * (ea - eb) / (a - b);

			if ((d1 * d1 - ga * gb) >= 0) {
				d2 = Math.sqrt(d1 * d1 - ga * gb);
				// Gpr.debug("gb: " + gb + "\tga: " + ga + "\td2: " + d2);
				alphaNew = b - (b - a) * (gb + d2 - d1) / (gb - ga + 2 * d2);
			} else alphaNew = a; // Forcing bisection

			interval = Math.abs(b - a);

			if ((Math.abs(a - alphaNew) < (interval * interSafeGuard)) || (Math.abs(b - alphaNew) < (interval * interSafeGuard))) alphaNew = (a + b) / 2; //If the cubic minimum is to close to the edges, the bisection method is choosen

			// Continue with zoom - calculating the properties of the new found step length
			for (i = 0; i < n; i++)
				x[i] = x0[i] + alphaNew * pk[i];

			energy.needsUpdate();
			eNew = energy.evaluate();
			if (debug) Gpr.debug("numAlphaEvaluations: " + numAlphaEvaluations + "\talphaNew: " + alphaNew + "\tx: " + Gpr.toString(x) + "\tx0: " + Gpr.toString(x0) + "\teNew: " + eNew);

			gradNew = 0; // Calculating the gradient at a=alphaNew
			double grad[] = energy.getGradient();
			for (i = 0; i < n; i++)
				gradNew += pk[i] * grad[i];

			if ((eNew > (e0 + c1 * alphaNew * grad0)) || (eNew >= eLow)) {
				alphaHi = alphaNew;
				eHi = eNew;
				gradHi = gradNew;
			} else {
				// This is the second Wolfe condition: p_k * grad[ f(x_k + alpha * p_k) ] >= c2 p_k * grad[ f(x_k) ]
				if (debug) Gpr.debug("Second Wolfe condition II:" //
						+ "\n\tp_k * grad[ f(x_k + alpha * p_k) ] >= c2 p_k * grad[ f(x_k) ] : " + (gradNew >= (c2 * grad0)) //
						+ "\n\tp_k * grad[ f(x_k + alpha * p_k) ] : " + gradNew //
						+ "\n\tc2 * p_k * grad[ f(x_k) ]          : " + (c2 * grad0) //
				);

				if (gradNew >= (c2 * grad0)) stop = true;
				else {
					if ((gradNew * (alphaHi - alphaLow)) >= 0) {
						alphaHi = alphaLow;
						eHi = eLow;
						gradHi = gradLow;
					}

					alphaLow = alphaNew;
					eLow = eNew;
					gradLow = gradNew;
				}

			}
		}

		alphaFinal = alphaNew;
	}
}
