 package meshi.PDB;
import meshi.util.filters.*;
public abstract class PdbLineFilter implements Filter {
    public boolean accept(Object obj) {
	return acceptPdbLine((PdbLine) obj);
    }
    public abstract boolean acceptPdbLine(PdbLine line);
}
