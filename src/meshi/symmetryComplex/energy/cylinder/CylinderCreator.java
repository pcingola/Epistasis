package meshi.symmetryComplex.energy.cylinder;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.Key;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;

/**
 * This energy function penalizes atoms which are outside a given cylinder.
 * More precisely, a "cylinder" is the space between two coaxial cylinders
 * with equal height, centered at the origin.
 * 
 * @version 0.1
 * @author Oren Wolfshtat
 */
public class CylinderCreator extends EnergyCreator implements KeyWords {
	protected double height = 0.0;
	protected double innerR = 0.0;
	protected double outerR = 0.0;
	protected Filter filter;

	public CylinderCreator() {
		super(CYLINDER_ENERGY);
	}

	public CylinderCreator(double weight) {
		super(weight);
	}

	public CylinderCreator(
			double outerR, double innerR, double height, Filter filter) {

		super(CYLINDER_ENERGY);

		if (innerR >= outerR)
			throw new RuntimeException(
			"CylinderCreator: innerR must be smaller than outerR");

		this.outerR = outerR;
		this.innerR = innerR;
		this.height = height;
		this.filter = filter;
	}

    public CylinderCreator(
			double outerR, double innerR, double height) {

		super(CYLINDER_ENERGY);

		if (innerR >= outerR)
			throw new RuntimeException(
			"CylinderCreator: innerR must be smaller than outerR");

		this.outerR = outerR;
		this.innerR = innerR;
		this.height = height;
	}

    public CylinderCreator(CommandList commands, String limits, Filter filter) {
		super(CYLINDER_ENERGY);
		Command limitsCommand = commands.firstWordFilter(CYLINDER_ENERGY).secondWord(limits);
		this.outerR = limitsCommand.thirdWordDouble();
		this.innerR = limitsCommand.fourthWordDouble();
		this.height = limitsCommand.fifthWordDouble();
		this.filter = filter;
		
		if (innerR >= outerR)
			throw new RuntimeException(
			"CylinderCreator: innerR must be smaller than outerR");
	}

	public AbstractEnergy createEnergyTerm(
			Protein protein, DistanceMatrix distanceMatrix,
			CommandList commands) {

		return new CylinderEnergy(
				protein.atoms(),
				distanceMatrix,
				new CylinderParametersList(height, innerR, outerR),
				weight(),
				filter);
	}

}