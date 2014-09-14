package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 09/06/2009
 * Time: 10:08:17
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZPropensityCreator extends EnergyCreator implements KeyWords {

    private static CooperativeZPropensityParameters parameters = null;
    private CompositePropensityCreator compositePropensityCreator;

    public CooperativeZPropensityCreator(CompositePropensityCreator compositePropensityCreator) {
		super(COOPERATIVE_Z_PROPENSITY_ENERGY );
		this.compositePropensityCreator = compositePropensityCreator;
	}


	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands ) {
		/* load parameters if this is the first time the function is called */
		if( parameters == null )
			parameters = new CooperativeZPropensityParameters( commands );
		CompositePropensityEnergy compositePropensityEnergy = (CompositePropensityEnergy) compositePropensityCreator.term();
		if (compositePropensityEnergy == null) throw new RuntimeException("Apparently an attempt to insantiate CooperativeZPropensity "+
										    "before CompositePropensityEnergy");
		return term = new CooperativeZPropensityEnergy(compositePropensityEnergy , weight(), parameters );
	}

}


