package meshi.energy.simpleEnergyTerms.alphaAngle;
import meshi.energy.simpleEnergyTerms.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;
import java.util.*;
             

public class AlphaAngleParametersList extends ParametersList {

    public AlphaAngleParametersList(String parametersFileName) {
	    super(parametersFileName,false);  // non-sortable
	}

    public Parameters createParameters(String line) {
    	return new AlphaAngleParameters(line);	
    }
    
    public Parameters parameters(Object obj) {
	Angle angle = (Angle) obj;
	Residue residue = angle.atom2.residue();
	int resnum = angle.getAngleResNum();
	String name = angle.getAngleName();
    String resName = angle.getAngleResName();
	AlphaAngleParameters alphaAngleParameters;

	if (! name.equals("CA3")) return null; 
	if (resnum < 0) return null; 	

	for (int cc=0 ; cc<size() ; cc++) {
		alphaAngleParameters = (AlphaAngleParameters) get(cc);
	    if (resName.equals(alphaAngleParameters.aaLetter)) 
	        return alphaAngleParameters;
	}

	return null;
    }
}
