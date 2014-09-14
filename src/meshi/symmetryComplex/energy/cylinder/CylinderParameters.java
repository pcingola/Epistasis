package meshi.symmetryComplex.energy.cylinder;

import meshi.energy.Parameters;

/**
 * @version 0.1
 * @author Oren Wolfshtat
 */
public class CylinderParameters implements Parameters {

	protected double height = 0.0;
	protected double innerR = 0.0;
	protected double outerR = 0.0;

	public CylinderParameters(
			double height, double innerR, double outerR) {

		this.height = height;
		this.innerR = innerR;
		this.outerR = outerR;
	}
}
