package meshi.sequences;
import java.io.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.*;

public class AtomAlignmentColumn extends AlignmentColumn{
    public AtomAlignmentColumn(Atom atom0, String comment0, Atom atom1, String comment1) {
	super(new AtomAlignmentCell(atom0, comment0),
	      new AtomAlignmentCell(atom1, comment1));
    }
    public AtomAlignmentColumn(Atom atom0, Atom atom1) {
	this(atom0, "protein0", atom1, "protein1");
    }
    
    Atom atomAt(int index) {return (Atom) cell(index).object();}
}
    
