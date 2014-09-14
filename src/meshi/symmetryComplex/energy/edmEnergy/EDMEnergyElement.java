package meshi.symmetryComplex.energy.edmEnergy;

import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.energy.EnergyElement;

public class EDMEnergyElement extends EnergyElement {

	private Atom atom;
	private EDMEnergyParameters ep;
	private double weight;
	
	public EDMEnergyElement( Atom atom, EDMEnergyParameters ep, double weight ) {
		this.atom = atom;
		this.ep = ep;
		this.weight = weight;
		
		setAtoms();
		updateFrozen();
	}
	
	public double evaluate() {
		if( frozen() ) return 0.0;
		
		double energy = ep.evaluate( 0, atom.x(), atom.y(), atom.z() ) * weight;
		double[] derivs = {
				ep.evaluate( 1, atom.x(), atom.y(), atom.z() ) * weight,
				ep.evaluate( 2, atom.x(), atom.y(), atom.z() ) * weight,
				ep.evaluate( 3, atom.x(), atom.y(), atom.z() ) * weight };
		
		atom.addToFx( derivs[0] );
		atom.addToFy( derivs[1] );
		atom.addToFz( derivs[2] );
		
		return energy;
	}

	protected void setAtoms() {
		atoms = new AtomList();
		atoms.add( atom );
	}

}
