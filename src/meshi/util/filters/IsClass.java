package meshi.util.filters;
public class IsClass implements Filter {
    private Class cls;
    public IsClass(Class cls) {
	this.cls = cls;
    }

    public boolean accept(Object obj) {
	return obj.getClass().equals(cls);
    }
} 
