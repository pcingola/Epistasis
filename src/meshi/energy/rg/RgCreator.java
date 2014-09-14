package meshi.energy.rg;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 *An implicit solvation energy term for all-atom models, modeling a 4.0 angs solvation shell around
 *each atom.
 **/

public class RgCreator extends EnergyCreator implements KeyWords {

	public RgCreator() {
		super(1.0);
	}

	public RgCreator(double weight) {
		super(weight);
	}

	@Override
	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
		return new RgEnergy(protein.atoms().filter(new AtomList.BackboneFilter()), distanceMatrix, weight());
	}

}
