package meshi.energy.simpleEnergyTerms.plane;
import java.util.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.util.string.*;
import meshi.geometry.*;
import meshi.util.*;

public class PlaneCreator extends SimpleEnergyCreator  implements KeyWords {
    public PlaneCreator(double weight) {
	super(weight);
    }
    public PlaneCreator() {
	super(PLANE_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new PlaneParametersList(parametersDirectory(commands)+
						     "/"+PLANE_PARAMETERS);
	
	TorsionList torsionList = TorsionList.createQuickAndDirtyTorsionList(protein, distanceMatrix);
	TorsionList relevantTorsionList = torsionList.filter(new HaveParametersFilter(parametersList));
    // Filtering out a problematic plane torsion in ARG (CD-NE-CZ-NH2)
    TorsionList finalTorsionList = new TorsionList();
    for (Torsion tor:relevantTorsionList)
       if (!(tor.atom1.name().equals("CD") && tor.atom2.name().equals("NE") &&
           tor.atom3.name().equals("CZ") && tor.atom4.name().equals("NH2")))
              finalTorsionList.add(tor);
	return new PlaneEnergy(finalTorsionList, distanceMatrix, (PlaneParametersList) parametersList, weight());
    }
}
	    
