package meshi.util.filters;
import meshi.util.*;
public class Differs implements Filter {
    private Object pivot;
    public Differs(Object obj) {
	pivot = obj;
    }
    public boolean accept(Object obj) {
	return (! pivot.equals(obj));
    }
}
