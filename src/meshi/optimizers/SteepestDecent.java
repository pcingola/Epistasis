package meshi.optimizers;

import meshi.optimizers.exceptions.LineSearchException;
import meshi.optimizers.exceptions.OptimizerException;

/**
 * This class implements a simple steepest descent minimizer, using a simple back tracking line search.
 * The search direction is given by the gradient of the energy function, at the minimizer position.
 * The step length first tried is the previous iteration step length, multiplied by some expansion factor.
 * If it leads to energy reduction it is chosen, otherwise it is shortened by some factor and tried again. This repeats
 * until energy reduction is achieved.
 *
 *To run this minimizer:
 *a) Instantiate this class with the desired minimization parameters.
 *b) Put the initial coordinates in the 'coordinates' variable at the 'energy' class.
 *c) Activate SteepestDecent.run().
 *d) Check for thrown errors to see if the minimization succeeded.
 *e) The minimized position is in the 'coordinates' variable at the 'energy' class.
 *
 * Parameters in the full constructor:
 * ----------------------------------
 * - energy - Energy object, where the energy function is.
 * - tolerance - Minimization stops when the magnitude of the maximal gradient component drops below tolerance.
 * - maxIteration - The maximal number of iteration steps allowed
 * - reoprtEvery - The frequency of the minimization reports.
 * - initialStepLength - parameter of the line search. The first step length to be tried after the calculation of the
 *                       first gradient. This parameter should normally be 1 unless very large gradients (such as clashing
 *                       of VDW atoms) are expected in the first steps. In that case it should be a much smaller number.
 * - stepSizeReduction - parameter of the line search. The step length is multiplied by this factor if no reduction
 *                       in energy is achieved.
 * - stepSizeExpansion - parameter of the line search. The first step length tried is the step length from previous
 *                       line search multiplied by this factor. (Note that non-positive values to this paramater cause
 *                       special options to be called (see the SimpleStepLength class help).
 **/

public class SteepestDecent extends Minimizer {

	public static final int DEFAULT_MAX_STEP = 1000;
	private SimpleStepLength lineSearch;

	private double lastStepLength = 1;
	private double initialStepLength;
	private double stepSizeReduction;
	private double stepSizeExpansion;

	// Default values constructor
	public SteepestDecent(Energy energy) {
		this(energy, DEFAULT_MAX_STEP, SimpleStepLength.DEFAULT_INITIAL_STEP_LENGTH, SimpleStepLength.DEFAULT_STEP_SIZE_REDUCTION, SimpleStepLength.DEFAULT_STEP_SIZE_EXPANTION);
	}

	//Full constructor
	public SteepestDecent(Energy energy, int maxIteration, double initialStepLength, double stepSizeReduction, double stepSizeExpansion) {
		super(energy);
		this.initialStepLength = initialStepLength;
		this.stepSizeReduction = stepSizeReduction;
		this.stepSizeExpansion = stepSizeExpansion;
		optimizerTerminator.setMaxSteps(maxIteration);
	}

	@Override
	protected void init() throws OptimizerException {
		lineSearch = new SimpleStepLength(energy(), initialStepLength, stepSizeReduction, stepSizeExpansion);
		lineSearch.setDebug(debug);
		energy().evaluate();
	}

	@Override
	protected void kickStart() throws OptimizerException {
	}

	public double lastStepLength() {
		return lastStepLength;
	}

	@Override
	protected boolean minimizationStep() throws OptimizerException {
		energy.evaluate();

		try {
			lastStepLength = lineSearch.findStepLength();
			if (lastStepLength <= 0) return true; // Could not improve after line search?
		} catch (LineSearchException lsException) {
			if (lsException.code == LineSearchException.NOT_A_DESCENT_DIRECTION) throw new OptimizerException("\n\nProblem in SteepestDecent. Direction is not a descent direction.\nThis problem is caused by incorrect differentiation of the energy function.\n");
			else throw new RuntimeException("Unknown LineSearchException " + lsException);
		}

		return true;
	}
}
