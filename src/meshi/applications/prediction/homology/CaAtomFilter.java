package meshi.applications.prediction.homology;

import meshi.util.filters.Filter;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;

public class CaAtomFilter implements Filter {
	public boolean accept(Object obj) {
	    return ((Atom) obj).name.equals("CA");
	}
}
