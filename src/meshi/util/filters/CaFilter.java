package meshi.util.filters;
import meshi.molecularElements.atoms.*;

public class CaFilter implements Filter {
    public static final CaFilter filter = new CaFilter();
    public boolean accept(Object obj) {
	if (((Atom) obj).name().equals("CA")) return true;
	return false;
    }
}
