package meshi.symmetryComplex.energy;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;

// TODO Remove abstract and constructor
public class SymmetryCreator extends EnergyCreator implements KeyWords{
	
	public SymmetryCreator() {
		super(SYMMETRY_ENERGY);
	}
	
	public SymmetryCreator(double weight) {
		super(weight);
	}
	
	public AbstractEnergy createEnergyTerm(
			Protein protein, DistanceMatrix distanceMatrix,
			CommandList commands) {
				
		if (! (protein instanceof SymmetricComplex))
			throw new RuntimeException
			("SymmetryCreator.createEnergyTerm -- first parameter must "+
			 "be an instance of SymmetricProtein.");

        
        ((SymmetricComplex)protein).resetNumberOfUpdates();
		return new SymmetryEnergy((SymmetricComplex)protein, distanceMatrix, null, weight());
		
//		return null;
	}
}