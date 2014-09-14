package meshi.energy.simpleEnergyTerms.plane;
import meshi.util.*;
import meshi.parameters.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import java.util.*;
             
public class PlaneParametersList extends ParametersList {
        public PlaneParametersList(String parametersFileName) {
	    super(parametersFileName,true);
	}

     public Parameters createParameters(String line) {
	return new PlaneParameters(line);
    }
    public Parameters parameters(Object obj) {
	Torsion torsion = (Torsion) obj;
	Parameters key = new PlaneParameters(torsion.atom1.type(), torsion.atom2.type(),
					     torsion.atom3.type(), torsion.atom4.type());
	return getParameters(key);
    }
}
