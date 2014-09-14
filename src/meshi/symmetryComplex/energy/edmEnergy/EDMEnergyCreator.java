package meshi.symmetryComplex.energy.edmEnergy;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;

import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;

/** The EDM energy class loads a splined polynomial for an electron
 * density map and attempts to fit the protein to that map.
 * 
 * @author El-ad David Amir
 */
public class EDMEnergyCreator	extends EnergyCreator
	implements KeyWords, MeshiPotential{
	
	protected Filter filter = null;

	public EDMEnergyCreator() {
		super( EDM_ENERGY );
	}
	
	public EDMEnergyCreator( double weight ) {
		super( weight );
	}
	
	public EDMEnergyCreator( Filter filter ) {
		super( EDM_ENERGY );
		this.filter = filter;
	}
	
	public EDMEnergyCreator( double weight, Filter filter ) {
		super( weight );
		this.filter = filter;
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		/* retrieve splined polynomial file name */
		String edmDataFileName = 
			commands.firstWord( EDM_ENERGY_FILE_NAME ).secondWord();

		/* create energy */
		return new EDMEnergy( protein,
				new EDMEnergyParametersList( edmDataFileName ),
				weight(),
				"EDM",
				filter);
	}
}


