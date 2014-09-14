package meshi.energy.simpleEnergyTerms.alphaTorsion;
import meshi.energy.simpleEnergyTerms.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;
import meshi.parameters.*;
import java.util.*;
             
public class AlphaTorsionParametersList extends ParametersList implements MeshiPotential {

    public AlphaTorsionParametersList(String parametersFileName) {
	    super(parametersFileName,false);  // non-sortable
	}

    public Parameters createParameters(String line) {
    	return new AlphaTorsionParameters(line);	
    }
    
    public Parameters parameters(Object obj) {
	Torsion torsion = (Torsion) obj;
	Residue residue = torsion.atom2.residue();
	int resnum = torsion.getTorsionResNum();
	String name = torsion.getTorsionName();
    String resName = torsion.getTorsionResName();
	AlphaTorsionParameters alphaTorsionParameters;

	if (! name.equals("ALPHA")) return null; 
	if (resnum < 0) return null; 	
	if (!((residue.secondaryStructure().equals(HELIX)) | 
	    (residue.secondaryStructure().equals(SHEET)))) return null;
 
	for (int cc=0 ; cc<size() ; cc++) {
		alphaTorsionParameters = (AlphaTorsionParameters) get(cc);
	    if (resName.equals(alphaTorsionParameters.aaLetter)) 
	        return alphaTorsionParameters;
	}
    
	return null;
    }
}
