package meshi.molecularElements.atoms;
import meshi.molecularElements.*;
import java.util.*;
import meshi.util.*;
import meshi.util.filters.*;

public class AtomPairList extends ArrayList<AtomPair> {
    /**
     * An empty AtomPair list
     **/
    public AtomPairList() {
	super();
    }
    public AtomPairList(ChainList chainList){
	this();
	Iterator bonds;
	AtomPair bond;
	ResidueList dummies = new ResidueList();
	for (Iterator chains = chainList.iterator(); chains.hasNext();){
	    Chain chain = (Chain) chains.next();
            Residue prevResidue = null;
 	    boolean nonDummyFound = false;
	    for (Residue residue:chain) {
		if (residue.dummy() &  nonDummyFound) dummies.add(residue); // this is a chain break
		if (! residue.dummy()) {
		    nonDummyFound = true;
		    add(residue.bonds());
		    if ((prevResidue != null) && (!prevResidue.dummy())){// we are not on the first residue.
			if (prevResidue.ID().number() != residue.ID().number() -1)
			    throw new RuntimeException("Weird consecutive residues "+prevResidue+" and "+residue);
		    
			if (prevResidue.ID().chain() != residue.ID().chain())
			    throw new RuntimeException("Weird consecutive residues "+prevResidue+" of chain "+
						       prevResidue.ID().chain()+
						       " and "+residue+" of chain "+residue.ID().chain());
			if ((prevResidue.nextAtom() == residue.tail()) &&
			    (residue.prevAtom() == prevResidue.head())) 
			    add(new AtomPair(prevResidue.head(),residue.tail()));
		        else {
				if (prevResidue.nextAtom() != null)   // Previous residue is already bound.
			    	throw new RuntimeException("Binding error 1 - An attempt to bind "+prevResidue+" to "+residue+
						       " while it is already bound to "+prevResidue.nextAtom());
				if (residue.prevAtom() != null)   // Current aresidue is already bound.
			    	throw new RuntimeException("Binding error 2 - An attempt to bind "+residue+" to "+prevResidue+
						       " while it is already bound to "+residue.prevAtom());
				if (prevResidue.head() == null) 
			    	throw new RuntimeException("Binding error 3 - An attempt to bind "+residue+" to "+prevResidue+
						       " while the head atom of "+prevResidue+" is null");
				if (residue.tail() == null) 
			    	throw new RuntimeException("Binding error 4 - An attempt to bind "+residue+
						       " with null tail to "+prevResidue);
				add(prevResidue.head().bond(residue.tail()));
				prevResidue.setNextAtom(residue.tail());
				residue.setPrevAtom(prevResidue.head());
			}
		    }
		}
		prevResidue = residue;
	    }
	}
   }
	    
    
    public AtomList atomList() {
	AtomList list = new AtomList();
	Iterator iter = iterator();
	Atom atom1, atom2;
	for (AtomPair atomPair:this){
	    atom1 = atomPair.atom1();
	    atom2 = atomPair.atom2();
	    if (! list.contains(atom1)) {
		list.add(atom1);
	    }
	    if (! list.contains(atom2)) {
		list.add(atom2);
	    }
	}
	return list;
    }
    public static class IsAtomPair implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof AtomPair); 
	}
    }

    public void printLaconic(String prompt) {
	for (AtomPair atomPair:this){
	    System.out.println(atomPair.laconic(prompt));
	}	
    }

    public void sort() {
	AtomPair[] array = toArray(new AtomPair[size()]);
	Arrays.sort(array);
	clear();
	for(AtomPair ap:array)
	    add(ap);
    }
    
    public void add(AtomPairList list) {
	for(AtomPair ap:list) add(ap);
    }

    public AtomPairList filter(Filter filter) {
	AtomPairList out = new AtomPairList();
	for (AtomPair ap:this) {
	    if (filter.accept(ap)) out.add(ap);
	}
	return out;
    }

}
    