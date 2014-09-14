package meshi.util.filters;
public class KolDichfin implements Filter {
	public boolean accept(Object obj) { 
		if (obj == null) return false;
		return true;
	}
}

