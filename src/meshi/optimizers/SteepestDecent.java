package meshi.optimizers;

import meshi.energy.Energy;
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
 *                       first gradient. This parameter should normally be 1 unless very large gradients (such as clashhing
 *                       of VDW atoms) are expected in the first steps. In that case it should be a much smaller number.
 * - stepSizeReduction - parameter of the line search. The step length is multiplied by this factor if no reduction
 *                       in energy is achieved.
 * - stepSizeExpansion - parameter of the line search. The first step length tried is the step length from previous
 *                       line search multiplied by this factor. (Note that non-positive values to this paramater cause
 *                       special options to be called (see the SimpleStepLength class help).
 **/

public class SteepestDecent extends Minimizer {

	private SimpleStepLength lineSearch;
	private double[][] coordinates;
	private double[][] bufferCoordinates;
	private double lastStepLength = 1;
	private static final double DEFAULT_TOLERANCE = 0.00001;
	private static final int DEFAULT_MAX_ITERATION = 100000;
	private static final int DEFAULT_REPORT_EVERY = 1;
	private static final double DEFAULT_INITIAL_STEP_LENGTH = 0.00000001;
	private static final double DEFAULT_STEP_SIZE_REDUCTION = 0.5;
	private static final double DEFAULT_STEP_SIZE_EXPENTION = 1.1;
	private double initialStepLength;
	private double stepSizeReduction;
	private double stepSizeExpansion;

	// Default values constructor
	public SteepestDecent(Energy energy) {
		this(energy, DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATION, DEFAULT_REPORT_EVERY, DEFAULT_INITIAL_STEP_LENGTH, DEFAULT_STEP_SIZE_REDUCTION, DEFAULT_STEP_SIZE_EXPENTION);
	}

	// Values for the line search are taken as default in this constructor
	public SteepestDecent(Energy energy, double tolerance, int maxSteps, int reportEvery) {
		this(energy, tolerance, maxSteps, reportEvery, DEFAULT_INITIAL_STEP_LENGTH, DEFAULT_STEP_SIZE_REDUCTION, DEFAULT_STEP_SIZE_EXPENTION);
	}

	//Full constructor
	public SteepestDecent(Energy energy, double tolerance, int maxSteps, int reportEvery, double initialStepLength, double stepSizeReduction, double stepSizeExpansion) {
		super(energy, maxSteps, reportEvery, tolerance);
		this.initialStepLength = initialStepLength;
		this.stepSizeReduction = stepSizeReduction;
		this.stepSizeExpansion = stepSizeExpansion;
	}

	@Override
	protected void init() throws OptimizerException {
		coordinates = energy().coordinates();
		lineSearch = new SimpleStepLength(energy(), initialStepLength, stepSizeReduction, stepSizeExpansion);
		bufferCoordinates = new double[coordinates.length][2];
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
		for (int i = 0; i < coordinates.length; i++) {
			bufferCoordinates[i][0] = coordinates[i][0];
			bufferCoordinates[i][1] = coordinates[i][1];
		}
		try {
			lastStepLength = lineSearch.findStepLength(bufferCoordinates);
		} catch (LineSearchException lsException) {
			if (lsException.code == LineSearchException.NOT_A_DESCENT_DIRECTION) throw new OptimizerException("\n\n Problem in SteepestDecent." + " Direction is not a descent direction. \n" + "This problem is caused by incorrect differentiation of the " + "energy function.\n" + "gradient rms is ");
			else throw new RuntimeException("Unknown LineSearchException " + lsException);
		}
		return true;
	}

	@Override
	public String toString() {
		return ("SteepestDecent\n" + "\t maxSteps \t" + maxSteps + "\n" + "\t tolerance \t" + tolerance);
	}
}
