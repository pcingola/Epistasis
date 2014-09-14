package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran;

import meshi.energy.EnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.parameters.*;

public class RamachandranEnergyElement
	extends EnergyElement
	implements CompositeTorsionsDefinitions {

	public final ResidueTorsions residueTorsions;
	public final RamachandranParameters cpp;
	protected double weight;

	public RamachandranEnergyElement(
			ResidueTorsions residueTorsions,
			RamachandranParameters cpp,
			double weight ) {
		this.residueTorsions = residueTorsions;
		this.cpp = cpp;
		this.weight = weight;
		
		setAtoms();
		updateFrozen();
	}

	protected void setAtoms() {
		int[] interestingTorsions = {PHI,PSI};
		atoms = residueTorsions.getAtoms(interestingTorsions);		
	}
	
	public double evaluate() {
		/* verify energy element is not frozen */
		if( frozen() ) return 0.0;
		
		/* calcualte energy and derivative */
		double energy       = cpp.evaluate( 0, residueTorsions );
		double phi_deriv	= cpp.evaluate( PHI, residueTorsions );
		double psi_deriv	= cpp.evaluate( PSI, residueTorsions );
		
		/* apply weight */
		energy *= weight;
		phi_deriv *= weight;
		psi_deriv *= weight;

		/* apply force to torsions */
		residueTorsions.applyForce( PHI, -phi_deriv );
		residueTorsions.applyForce( PSI, -psi_deriv );

		return energy;		
	}	
}
