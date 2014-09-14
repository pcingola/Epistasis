package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;

import meshi.energy.*;
import meshi.molecularElements.*;
import meshi.geometry.*;
import meshi.util.*;

public class CooperativeRamachandranCreator extends EnergyCreator implements KeyWords {
    
    private static CooperativeRamachandranParameters parameters = null;
    private RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator;
	public CooperativeRamachandranCreator(RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator) {
		super(COOPERATIVE_RAMACHANDRAN_ENERGY );
		this.ramachandranSidechainEnergyCreator = ramachandranSidechainEnergyCreator;
	}


	public AbstractEnergy createEnergyTerm( Protein protein, DistanceMatrix distanceMatrix, CommandList commands ) {
		/* load parameters if this is the first time the function is called */
		if( parameters == null )
			parameters = new CooperativeRamachandranParameters( commands );
		RamachandranSidechainEnergy ramachandranSidechainEnergy = (RamachandranSidechainEnergy) ramachandranSidechainEnergyCreator.term();
		if (ramachandranSidechainEnergy == null) throw new RuntimeException("Apparently an attempt to insantiate CooperativeRamachandran "+
										    "before RamachandranSidechainEnergy");
		return term = new CooperativeRamachandran(ramachandranSidechainEnergy , weight(), parameters );
	}

}