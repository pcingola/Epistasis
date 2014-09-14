package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

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
 * Time: 13:21:12
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZRamachandranCreator extends EnergyCreator implements KeyWords {

    private static CooperativeZRamachandranParameters parameters = null;
    private RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator;
	public CooperativeZRamachandranCreator(RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator) {
		super(COOPERATIVE_Z_RAMACHANDRAN_ENERGY );
		this.ramachandranSidechainEnergyCreator = ramachandranSidechainEnergyCreator;
	}


	public AbstractEnergy createEnergyTerm( Protein protein, DistanceMatrix distanceMatrix, CommandList commands ) {
		/* load parameters if this is the first time the function is called */
		if( parameters == null )
			parameters = new CooperativeZRamachandranParameters( commands );
		RamachandranSidechainEnergy ramachandranSidechainEnergy = (RamachandranSidechainEnergy) ramachandranSidechainEnergyCreator.term();
		if (ramachandranSidechainEnergy == null) throw new RuntimeException("Apparently an attempt to insantiate CooperativeZRamachandran "+
										    "before RamachandranSidechainEnergy");
		return term = new CooperativeZRamachandran(ramachandranSidechainEnergy , weight(), parameters );
	}

}

