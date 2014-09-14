package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;
import meshi.util.*;

public class RamachandranSidechainEnergyCreator
	extends EnergyCreator
	implements KeyWords, MeshiPotential {

	public RamachandranSidechainEnergyCreator() {
		super( RAMACHANDRAN_SIDECHAIN_ENERGY );
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		/* retrieve parameters */
		String rsplFileName = 
			parametersDirectory(commands)+"/"+COMPOSITE_TORSIONS_PARAMETERS;
		RamachandranSidechainParametersList rspl = 
			new RamachandranSidechainParametersList(rsplFileName);
		
		/* create residue torsions list for protein */
		ResidueTorsionsList rtl = new ResidueTorsionsList(
				protein, distanceMatrix );
		
		/* return energy */
		return term = new RamachandranSidechainEnergy(rtl, distanceMatrix, rspl, weight(), "EnResidue" );
	}

}
