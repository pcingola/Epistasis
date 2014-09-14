package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;
import meshi.parameters.*;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/** The Ramachandran energy calculates the probability of each amino acid to be in a specific position 
 *  at the (phi,psi) torsion space. 
 */
public class RamachandranCreator extends EnergyCreator implements KeyWords, MeshiPotential{
    private static RamachandranParametersList parametersList = null;
	public RamachandranCreator(double weight) {
	super(weight);
    }
    
    public RamachandranCreator() {
		super( 1.0 );
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		/* retrieve parameters */
		if (parametersList == null) {
			String cpplFileName = parametersDirectory(commands)+"/"+COMPOSITE_PROPENSITY_2D_PARAMETERS;
			parametersList  = new RamachandranParametersList(cpplFileName);
		}
		
		/* create residue torsions list for protein */
		ResidueTorsionsList rtl = (new ResidueTorsionsList(protein, distanceMatrix)).filterPhiPsiResidues();
		
		/* return energy */
		return new RamachandranEnergy(rtl, distanceMatrix, 
		(RamachandranParametersList) parametersList, weight(), "Ramach" );
	}	
}
