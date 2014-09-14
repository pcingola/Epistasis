package meshi.util;
public class MeshiException extends RuntimeException {
    public MeshiException(String errorMessage) {
	super();
	System.err.print("\n"+errorMessage+"\n");
	if (MeshiProgram.debug()) 
	   MeshiProgram.printGlobalTable();
	printStackTrace();
    }
}
