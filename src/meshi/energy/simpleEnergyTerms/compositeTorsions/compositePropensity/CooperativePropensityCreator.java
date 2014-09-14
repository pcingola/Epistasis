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
 * Date: 09/02/2009
 * Time: 10:45:41
 * To change this template use File | Settings | File Templates.
 */
public class CooperativePropensityCreator extends EnergyCreator implements KeyWords {

    private static CooperativePropensityParameters parameters = null;
    private CompositePropensityCreator compositePropensityCreator;

    public CooperativePropensityCreator(CompositePropensityCreator compositePropensityCreator) {
		super(COOPERATIVE_PROPENSITY_ENERGY );
		this.compositePropensityCreator = compositePropensityCreator;
	}


	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands ) {
		/* load parameters if this is the first time the function is called */
		if( parameters == null )
			parameters = new CooperativePropensityParameters( commands );
		CompositePropensityEnergy compositePropensityEnergy = (CompositePropensityEnergy) compositePropensityCreator.term();
		if (compositePropensityEnergy == null) throw new RuntimeException("Apparently an attempt to insantiate CooperativePropensity "+
										    "before CompositePropensityEnergy");
		return term = new CooperativePropensityEnergy(compositePropensityEnergy , weight(), parameters );
	}

}
