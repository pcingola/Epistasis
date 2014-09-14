package meshi.energy.simpleEnergyTerms.angle;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.string.*;
import meshi.util.*;

public class AngleCreator extends SimpleEnergyCreator  implements KeyWords {
    public AngleCreator() {
	super(ANGLE_ENERGY);
    }

    public AngleCreator(double weight) {
	super(weight);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
 	AtomPairList bondList = (AtomPairList) protein.bonds().filter(new GoodBonds());
	AngleList angleList = new QuickAndDirtyAngleList(bondList, distanceMatrix);
	if (parametersList== null)
	    parametersList = new AngleParametersList(parametersDirectory(commands)+
						     "/"+ANGLE_PARAMETERS);
	return new AngleEnergy(angleList, distanceMatrix, (AngleParametersList) parametersList, weight());
    }
}
