package meshi.energy.simpleEnergyTerms.bond;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.util.string.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.*;
/**
 * A factory for BondEnergy objects.
 **/
public class BondCreator extends SimpleEnergyCreator  implements KeyWords {
    public BondCreator() {
	super(BOND_ENERGY);
    }

    public BondCreator(double weight) {
	super(weight);
    }

    /**
     *<pre> 
     * hides all the hard work needed to generate a BondEnergy object. 
     * a) Extract the bonds from the protein.
     * b) Finds and reads the parameters file.
     **/ 
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
 	AtomPairList bondList = (AtomPairList) protein.bonds().filter(new GoodBonds());
	if (parametersList == null)
	    parametersList = new BondParametersList(parametersDirectory(commands)+
						    "/"+BOND_PARAMETERS);
	term = new BondEnergy(bondList, distanceMatrix, (BondParametersList) parametersList, weight());
	return term;
    }    
}
