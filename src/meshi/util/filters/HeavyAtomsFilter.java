package meshi.util.filters;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

public class HeavyAtomsFilter implements Filter {
    public boolean accept(Object obj) {
	Atom atom = (Atom) obj;
	if (atom.type().isHydrogen()) return false;
	return true;
    }
}
