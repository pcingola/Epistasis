package meshi.energy.simpleEnergyTerms.outOfPlane;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*; 
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.util.string.*;
import meshi.geometry.*;
import meshi.util.*;
import meshi.util.filters.*;

           
public class OutOfPlaneCreator extends SimpleEnergyCreator  implements KeyWords {
    public OutOfPlaneCreator(double weight) {
	super(weight);
    }
    public OutOfPlaneCreator() {
	super(OUT_OFPLANE_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new OutOfPlaneParametersList(parametersDirectory(commands)+
							  "/"+OUT_OF_PLANE_PARAMETERS);
	TorsionList torsionList = TorsionList.createQuickAndDirtyTorsionList(protein, distanceMatrix);
	    TorsionList OOPlist = (TorsionList) torsionList.filter(new TorsionList.FilterOOP());
	TorsionList relevantTorsionList = (TorsionList) OOPlist.filter(new HaveParametersFilter(parametersList));
	return new OutOfPlaneEnergy(distanceMatrix, relevantTorsionList, (OutOfPlaneParametersList) parametersList, weight());
    }
}
