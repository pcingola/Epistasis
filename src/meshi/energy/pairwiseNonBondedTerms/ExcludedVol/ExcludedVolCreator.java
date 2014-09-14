package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;

public class ExcludedVolCreator extends EnergyCreator  implements KeyWords {


    private int type=0; // allows us to run one of the EV scenario
    private double Rfac = 1.0;
    private Filter filter = null; 
    private static ExcludedVolParametersList parametersList = null;
  
    public ExcludedVolCreator(double weight , int type) {
  	super(weight);
  	this.type = type;
    }

    public ExcludedVolCreator(double weight) {
  	super(weight);
    }

    public ExcludedVolCreator() {
  	super(EXCLUDED_VOL);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new ExcludedVolParametersList(parametersDirectory(commands)+
							    "/"+EXCLUDED_VOL_PARAMETERS);
	return new ExcludedVol(distanceMatrix, parametersList, type,  weight(),Rfac,filter);
    }
}
