package meshi.sequences;
import meshi.util.filters.*;
import meshi.util.*;
import java.util.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

public class AtomAlignment extends ArrayList<AtomAlignmentColumn>{
    public AtomAlignment() {
	super();
    }
    public AtomAlignment(AtomList atomList0,AtomList atomList1) {
	this();
	if (atomList0.size() != atomList1.size()) throw new RuntimeException("List lengths must be equal");
	Iterator atoms1 = atomList1.iterator();
	for (Iterator atoms0 = atomList0.iterator(); atoms0.hasNext();)
	    add(new AtomAlignmentColumn((Atom) atoms0.next(), atomList0.comment(),
					(Atom) atoms1.next(), atomList1.comment()));
    }

    public AtomAlignment(ResidueAlignment residueAlignment) {
	    this(residueAlignment, new KolDichfin());
    }

    public AtomAlignment(ResidueAlignment residueAlignment, Filter filter) {
	super();
	if (residueAlignment.get(0).size() != 2) 
	    throw new RuntimeException("This constructor supports only pairwise residueAlignments");
	for(Iterator columns = residueAlignment.iterator(); columns.hasNext();) {
	    ResidueAlignmentColumn column = (ResidueAlignmentColumn) columns.next();
	    if (!column.hasGap()) {
		for (Iterator atoms0 = column.residue0().atoms().iterator(); atoms0.hasNext();) {
		    Atom atom0 = (Atom) atoms0.next();
		    if (filter.accept(atom0) & (!atom0.nowhere()))
			for (Iterator atoms1 = column.residue1().atoms().iterator(); atoms1.hasNext();) {
			    Atom atom1 = (Atom) atoms1.next();
			    if ((filter.accept(atom1)) & (atom0.name().equals(atom1.name())) & (!atom1.nowhere()))
				add(new AtomAlignmentColumn(atom0, atom1));
			}
		}
	    }
	}
    }
		

    private static class IsAlignmentColumn implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof AtomAlignmentColumn);
	}
    }
    
    public Rms rms() {
	if (hasGaps()) throw new RuntimeException("Cannot calculate RMS for AtomAlignment with gaps");
	return new Rms(this);
    }

    public Atom atomAt(int coulumnIndex, int rowIndex) {
	return (get(coulumnIndex)).atomAt(rowIndex);
    }
    
    public AtomList atomList(int row) {
	AtomList out = new AtomList();
	for (Iterator columns = iterator(); columns.hasNext();) {
 	    AtomAlignmentColumn  column = (AtomAlignmentColumn) columns.next();
 	    AtomAlignmentCell cell = (AtomAlignmentCell) column.cell(row);
	    out.add(cell.atom());
	}
	return out;
    }

    /**
     * Does the alignemnt include gaps? .
     **/
    public boolean hasGaps() { 
        for (Iterator iter = iterator(); iter.hasNext();) { 
            AlignmentColumn column = (AlignmentColumn) iter.next(); 
            if (column.hasGap()) return true; 
        } 
        return false; 
    } 
 }
