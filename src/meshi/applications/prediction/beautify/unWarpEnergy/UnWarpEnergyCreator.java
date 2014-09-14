package meshi.applications.prediction.beautify.unWarpEnergy;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;


public class UnWarpEnergyCreator extends EnergyCreator {

    public UnWarpEnergyCreator() {
        super(UN_WARP_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        term = new UnWarpEnergy(weight(), protein);
        return term;
    }
}
   
