package meshi.sequences;
import java.util.*;

public class AccesibilitySequence  extends Sequence {
    public static final AccesibilitySequenceCharFilter AccesibilityCharFilter =  new AccesibilitySequenceCharFilter();
    public AccesibilitySequence(String sequence, String comment) {
	super(sequence, comment, AccesibilityCharFilter);
    }
    
        private static class AccesibilitySequenceCharFilter extends SequenceCharFilter {
	public boolean accept(Object obj) {
	    Character c = ((Character) obj).charValue();
	    if ("AB".indexOf(c) >= 0) return true;
	    return false;
	}
    }
}
