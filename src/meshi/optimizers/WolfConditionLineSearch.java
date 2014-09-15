package meshi.optimizers;

import meshi.optimizers.exceptions.LineSearchException;

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

public class WolfConditionLineSearch extends LineSearch {

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
	private double alphaI, alphaI1, alphaFinal;
	private int n; // Number of variables
	private double e0, eI, eI1; // Energy values along Pk
	private double grad0, gradI, gradI1; // Gradients along Pk
	private int numAlphaEvaluations, stop;
	private int i;

	public WolfConditionLineSearch(Energy energy) {
		this(energy, DEFAULT_C1, DEFAULT_C2, DEFAULT_EXTENDED_ALPHA_FACTOR, DEFAULT_MAX_NUM_EVALUATIONS);
	}

	public WolfConditionLineSearch(Energy energy, double c1, double c2, double extendAlphaFactor, int maxNumEvaluations) {
		super(energy);
		this.c1 = c1;
		this.c2 = c2;
		this.maxNumEvaluations = maxNumEvaluations;
		this.extendAlphaFactor = extendAlphaFactor;
		interSafeGuard = extendAlphaFactor / interSafeGuardFactor;

		x = energy.getX(); // Local copy of coordinates
		n = x.length;
		xCopy = new double[n];

		reset();
	}

	@Override
	public double findStepLength() throws LineSearchException {
		energy.copyX(xCopy);
		return findStepLength(xCopy, energy.getGradient());
	}

	public double findStepLength(double x0[], double pk[]) throws LineSearchException {

		// Initializing the run
		energyNew = energy.getEnergy();

		// Calculating  grad_0 = - p_k^T * grad[ f(x_k) ]
		grad0 = 0;
		double grad[] = energy.getGradient();
		for (i = 0; i < n; i++)
			grad0 -= pk[i] * grad[i];

		// Sanity check: Is this in the descent direction?
		if (grad0 >= 0) throw new LineSearchException(LineSearchException.NOT_A_DESCENT_DIRECTION, "\n\nThe search direction is not a descent direction. \n" + "This problem might be caused by incorrect diffrentiation of the energy function,\n" + "or by numerical instabilities of the minimizing techniques" + "(such as not fullfilling the Wolf condtions in BFGS).\n");

		numAlphaEvaluations = -1;
		stop = 0;
		e0 = energyNew;
		eI = e0;
		alphaI = 0;

		// Choosing the initial alpha guess
		if (firstRun != 0) {
			firstRun = 0;
		} else {
			alphaI1 = 2.02 * (energyNew - energyOld) / grad0;
			if (alphaI1 > 1) alphaI1 = 1;
		}

		energyOld = energyNew; // For the next time line search is called

		// Bracketing the Wolf area
		while ((numAlphaEvaluations <= maxNumEvaluations) && (stop == 0)) {
			numAlphaEvaluations++;
			for (i = 0; i < n; i++)
				x[i] = x0[i] + alphaI1 * pk[i];

			eI1 = energy.evaluate();
			gradI1 = 0; // calculating the gradient at a=alphaI1
			for (i = 0; i < n; i++)
				gradI1 -= pk[i] * grad[i];

			if ((eI1 > (e0 + c1 * alphaI1 * grad0)) || ((eI1 >= eI) && (numAlphaEvaluations > 0))) zoom(x0, false); // Calling the regular zoom
			else {
				if (Math.abs(gradI1) <= (-c2 * grad0)) {
					alphaFinal = alphaI1;
					stop = 1;
				} else if (gradI1 >= 0) zoom(x0, true); // Inverse Zoom
			}

			alphaI = alphaI1;
			gradI = gradI1;
			eI = eI1;
			alphaI1 = alphaI1 * extendAlphaFactor;
		}

		if (numAlphaEvaluations > maxNumEvaluations) {
			energy.setX(x0); // Returning the coordinates to the original state
			energy.evaluate();
			throw new LineSearchException(LineSearchException.WOLF_CONDITION_NOT_MET, "\n\nWolf conditions not met. The line search did not converge, and exceeded the maximal number of step extensions allowed.\n");
		}

		return alphaFinal;
	}

	public void reset() {
		firstRun = 1;
		alphaI1 = 1.0;
	}

	public void reset(double resAlpha) {
		firstRun = 1;
		alphaI1 = resAlpha;
	}

	// The function Zoom finds a step length satisfing the Wolf conditions, given the bracketing of alphaI and alphaI1.
	// It was separated into a different function to make the code more readable.
	private void zoom(double[] x0, boolean inv) {
		double alphaHi, alphaLow, alphaNew = 0;
		double eHi, eLow, eNew = 0;
		double gradHi, gradLow, gradNew = 0;
		double a = 0, b = 0, ga = 0, gb = 0, ea, eb, interval, d1 = 0, d2 = 0; // Cubic interpolation variables

		eHi = eI1;
		gradHi = gradI1;
		alphaHi = alphaI1;
		eLow = eI;
		gradLow = gradI;
		alphaLow = alphaI;

		if (inv) {
			eHi = eI;
			gradHi = gradI;
			alphaHi = alphaI;
			eLow = eI1;
			gradLow = gradI1;
			alphaLow = alphaI1;
		}

		while ((numAlphaEvaluations <= maxNumEvaluations) && (stop == 0)) {

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
				alphaNew = b - (b - a) * (gb + d2 - d1) / (gb - ga + 2 * d2);
			} else alphaNew = a; // Forcing bisection

			interval = Math.abs(b - a);

			if ((Math.abs(a - alphaNew) < (interval * interSafeGuard)) || (Math.abs(b - alphaNew) < (interval * interSafeGuard))) alphaNew = (a + b) / 2; //If the cubic minimum is to close to the edges, the bisection method is choosen

			// Continue with zoom - calculating the properties of the new found step length
			for (i = 0; i < n; i++)
				x[i] = x0[i] + alphaNew * x0[i];

			eNew = energy.evaluate();

			gradNew = 0; // calculating the gradient at a=alphaNew
			for (i = 0; i < n; i++)
				gradNew -= x0[i] * x[i];

			if ((eNew > (e0 + c1 * alphaNew * grad0)) || (eNew >= eLow)) {
				alphaHi = alphaNew;
				eHi = eNew;
				gradHi = gradNew;
			} else {

				if (Math.abs(gradNew) <= (-c2 * grad0)) stop = 1;
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
