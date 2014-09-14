package meshi.parameters;
import meshi.util.*;

public enum SecondaryStructure implements MeshiAttribute{
    HELIX("H"),
    SHEET("E"),
    COIL("C"),
    HELIX_OR_COIL("h"),
    SHEET_OR_COIL("e"),
    ALL("A"),
    UNK("X");

    private final String nameOneLetter;
    private SecondaryStructure(String nameOneLetter) {
	this.nameOneLetter = nameOneLetter;
    }
    
    public static SecondaryStructure secondaryStructure(char c) {
	for (SecondaryStructure sec:SecondaryStructure.values()) 
	    if (sec.nameOneLetter.charAt(0)==c) return sec;
	throw new RuntimeException("Undefined secondary structure: "+c+"\n"+"Please take a look at meshi.parameters.SecondaryStructure.");
    }
    public static boolean isSecondaryStructure(char c) {
	for (SecondaryStructure sec:SecondaryStructure.values()) 
	    if ((sec.nameOneLetter.charAt(0)==c) &
		(sec!=ALL) &
		(sec!=UNK))return true;
	return false;
    }

    public int key() {
	return SECONDARY_STRUCTURE_ATTRIBUTE;
    }
}
