package meshi.molecularElements;
import meshi.molecularElements.atoms.*; 
import meshi.util.*;
import java.util.*;
import meshi.util.filters.Filter;
import meshi.parameters.*;
import meshi.sequences.*;

/**
 * A list of residues.
 * 
 **/
public class ResidueList extends ArrayList<Residue>{

    // --------------------------------- constructors ------------------------------
    
    public ResidueList() {
	super();
    }

    public ResidueList(ChainList chains) {
	super();  /*
        for (Iterator chainsIter = chains.iterator(); chainsIter.hasNext();) {
            Chain chain = (Chain) chainsIter.next();
            Iterator residues = chain.iterator();
            Residue residue = (Residue) residues.next();
            if (! residue.dummy())
            throw new RuntimeException("The first residue in each Chain(so as in Residues)  must be dummy!");
            add(residue);

            for (; residues.hasNext();) {
                residue = (Residue) residues.next();
                if (! residue.dummy()) add(residue);
            }
        }

        /*/
    for (Iterator chainsIter = chains.iterator(); chainsIter.hasNext();) {
	    Chain chain = (Chain) chainsIter.next();
	    for (Iterator residues = chain.iterator(); residues.hasNext();) {
		Residue residue = (Residue) residues.next();
		if (! residue.dummy()) add(residue);
	    }
	}     //*/

    }
    //------------------------------------------- methods -----------------------------------
    /**
     * Fetch a residue by its position in the list.
     **/
    
    /**
     * Fetches a residue by its residue identifier.
     **/
    public Residue residue(ResidueIdentifier residueID) {
	Residue residue;
	Iterator residueIterator = iterator();
	while (residueIterator.hasNext()) {
	    residue = (Residue) residueIterator.next();
	    if (residue.ID().equals(residueID))
		return residue;
	}
	return null;
    }



    /**
     * Extracts the residues accepted by the filter.
     **/
    public ResidueList filter(Filter residueFilter) {
	ResidueList newList = new ResidueList();
	for(Residue residue:this){
	    if (residueFilter.accept(residue))
		newList.add(residue);
	}
	return newList;
    }

  



	    

	public String toString() {
	String out = "";
	Iterator residues = iterator();
	boolean first = true;;
	while (residues.hasNext()) {
	    Residue residue = (Residue) residues.next();
	    if ((!first) || (!residue.dummy())) {
		if (residue.dummy()) out+=SequenceAlignmentCell.WILDCARD_CHAR;
		else out+=residue.type().nameOneLetter();
	    }
	    first = false;
	}
	return out;
    }

    public AtomList atoms() {return new AtomList(this);}


    public static class NonDummyFilter implements Filter {
	public boolean accept(Object obj) {
	    return !((Residue) obj).dummy();
	}
    }

    public void print() {
	for (Residue r:this)
	    System.out.println(r);
    }

    public void sort() {
	Residue[] temp = toArray(new Residue[size()]);
	clear();
	for (Residue r:temp)
	    add(r);
    }

}
