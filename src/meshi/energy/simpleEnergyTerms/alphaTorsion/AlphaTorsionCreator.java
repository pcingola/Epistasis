package meshi.energy.simpleEnergyTerms.alphaTorsion;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.util.string.*;
import meshi.geometry.*;
import meshi.util.*;

public class AlphaTorsionCreator extends SimpleEnergyCreator  implements KeyWords {
    public AlphaTorsionCreator(double weight) {
	super(weight);
    }
    public AlphaTorsionCreator() {
	super(ALPHA_TORSION_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
					  	
	if (parametersList== null) {                                 
	    parametersList = new AlphaTorsionParametersList(parametersDirectory(commands)+
						     "/"+ALPHA_TORSION_PARAMETERS);
	}
	QuickAndDirtyAngleList angleList  = QuickAndDirtyAngleList.getCaAnglesQuickAndDirty(protein,distanceMatrix).namedFilterQD();
	TorsionList torsionList = new TorsionList(angleList, distanceMatrix);
	TorsionList relevantTorsionList = (TorsionList)torsionList.filter(new HaveParametersFilter(parametersList)); 

	term = new AlphaTorsionEnergy(relevantTorsionList, distanceMatrix, (AlphaTorsionParametersList) parametersList, weight());
	return term;
    }
}
