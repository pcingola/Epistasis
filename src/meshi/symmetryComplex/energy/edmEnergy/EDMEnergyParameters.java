package meshi.symmetryComplex.energy.edmEnergy;

import meshi.energy.Parameters;

/** A dummy parameters class that always includes a splined polynomial. */
public class EDMEnergyParameters implements Parameters {

	private EDMSplinedPolynomial esp = null;
	
	public EDMEnergyParameters( EDMSplinedPolynomial esp ) {
		this.esp = esp;
	}

	/** Return value of polynomial at given points, derive given variable. */
	public double evaluate( int derivVar, double ... args ) {
		return esp.value( derivVar, args );
	}
		
}


