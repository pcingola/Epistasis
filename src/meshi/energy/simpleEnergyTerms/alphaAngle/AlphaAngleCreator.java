package meshi.energy.simpleEnergyTerms.alphaAngle;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.util.string.*;
import meshi.geometry.*;
import meshi.util.*;

public class AlphaAngleCreator extends SimpleEnergyCreator  implements KeyWords {
    public AlphaAngleCreator(double weight) {
	super(weight);
    }
    public AlphaAngleCreator() {
	super(ALPHA_ANGLE_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	QuickAndDirtyAngleList angleList  = QuickAndDirtyAngleList.getCaAnglesQuickAndDirty(protein,distanceMatrix).namedFilterQD();
	if (parametersList== null) {                                 
	    parametersList = new AlphaAngleParametersList(parametersDirectory(commands)+
						     "/"+ALPHA_ANGLE_PARAMETERS);
	 }
	term = new AlphaAngleEnergy(angleList, distanceMatrix, (AlphaAngleParametersList) parametersList, weight());
	return term;
    }
}
