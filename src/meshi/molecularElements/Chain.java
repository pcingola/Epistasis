package meshi.molecularElements;
import meshi.molecularElements.atoms.*; 
import meshi.util.*;
import java.util.*;
import meshi.util.filters.Filter;
import meshi.parameters.*;
import meshi.sequences.*;
import meshi.geometry.*;

/**
 * A list of residues of the same chain.
 * Note that some of the residues may be dummy. 
 * Specifically, the first residue may be dummy in order to be compatible with the 
 * biologists convention that the first residue is 1 (they grew up on FORTRAN). 
 * Residue positions in the list are equal to their residue number. Holes are filed by DummyResidue objects 
 **/
public class Chain extends ResidueList{
    /**
     * The first Non Dummy Residue (may but also may not be number 1).
     **/
    public final String name;
    public final String name() {return name;}
    public static final String GENERIC_CHAIN_NAME = " ";
    public final Protein protein;

    // --------------------------------- constructors ------------------------------
    /**
     * Empty list.
     **/    
    public Chain(String name, Protein protein) {
	super();
	this.name = name; 
	this.protein = protein;
    }

    /**
     * Constructs a Chain from a sequence (a string of one letter codes).
     * The residueCreator is responsible for the interpretation of the letters to actual residues and thus
     * determines the molecular model. 
     **/
    public Chain(Sequence sequence, ResidueCreator creator, Protein protein) {
	this(sequence, creator,GENERIC_CHAIN_NAME, protein);
    }

     /**
     * Constructs a Chain from a sequence (a string of one letter codes).
     * The residueCreator is responsible for the interpretation of the letters to actual residues and thus
     * determines the molecular model. 
     **/
     public Chain(Sequence sequence, ResidueCreator creator, String chainName, Protein protein) {
     this(chainName, protein);
     int i;
     Residue newResidue;
     SecondaryStructure ss;
     add(new Residue(new ResidueIdentifier(chainName,0))); //traditionally sequences start at 1
     if (sequence.size() < 1) throw new RuntimeException("Protein sequence length needs to be at least one");
     if (sequence.size() < 2)
         add(creator.create(getCaAtomList(sequence.getChar(0),
                          new Coordinates()),
                    new ResidueIdentifier(chainName,1),
                    ResidueMode.SINGLE));
     else {
         i = 0;
         while(((i < sequence.startsIn()-1) | (sequence.getChar(i) == 'X'))
           & (i < sequence.size())) {
         newResidue = new Residue(new ResidueIdentifier(chainName,i+1));
         add(newResidue);
         i++;
         }
         if (i == sequence.size()) throw new RuntimeException("Weired sequence");

         AtomList atomList = getCaAtomList(sequence.getChar(i-sequence.startsIn()+1),
                           new Coordinates());
         newResidue = creator.create(atomList,
                     new ResidueIdentifier(chainName,i+1),
                     ResidueMode.NTER);
         ss = sequence.getSs(i-sequence.startsIn()+1);
         newResidue.setSecondaryStructure(ss);
         add(newResidue);
         i++;
         for (; i < sequence.size()+sequence.startsIn()-2; i++) {
             if (sequence.getChar(i-sequence.startsIn()+1) == 'X'){
                 newResidue = new Residue(new ResidueIdentifier(chainName,i+1));
                 add(newResidue);
             }
             else {
                 newResidue = creator.create(getCaAtomList(sequence.getChar(i-sequence.startsIn()+1),
                                       new Coordinates()),
                                 new ResidueIdentifier(chainName,i+1),
                                 ResidueMode.NORMAL);
                 ss = sequence.getSs(i-sequence.startsIn()+1);
                 newResidue.setSecondaryStructure(ss);
                 add(newResidue);
             }
         }
         newResidue = creator.create(getCaAtomList(sequence.getChar(i-sequence.startsIn()+1),
                               new Coordinates()),
                     new ResidueIdentifier(chainName,i+1),
                     ResidueMode.CTER);
         ss = sequence.getSs(i-sequence.startsIn()+1);
         newResidue.setSecondaryStructure(ss);
         add(newResidue);
     }

     }
    
   /**
     * Constructs a Chain from a list of atoms.
     * The residueCreator allows the use of information that is not  stored in the atoms themselves. 
     * Note that the chain is going to havewill get a genergic name.
     **/
    public Chain(AtomList atomList, ResidueCreator creator, Protein protein) {
	this(atomList, creator, GENERIC_CHAIN_NAME, protein);
    }

   /**
     * Constructs a Chain from a list of atoms.
     * The residueCreator allows the use of information that is not  stored in the atoms themselves. 
     **/
    public Chain(AtomList atomList, ResidueCreator creator, String chainName, Protein protein) {
	this(chainName, protein);
	AtomList newAtomList = null;
	boolean first = true;
	add(new Residue(new ResidueIdentifier(chainName,0))); //traditionally sequences start at 1
	if (atomList.size() == 0) throw new RuntimeException(" No Atoms in AtomList "+atomList.comment());
	for (Atom atom:atomList){
	    while (size()< atom.residueNumber()) {
		if (newAtomList == null) add(new Residue(new ResidueIdentifier(chainName,size())));
		else {
		    if (first) {
				createAndAddResidue(creator, newAtomList, ResidueMode.NTER);
				first = false;
		    }
		    else {
				createAndAddResidue(creator, newAtomList, ResidueMode.NORMAL);
			}
			newAtomList = null;
		}
	    }
	    if (size() == atom.residueNumber()) {
		if (newAtomList == null) newAtomList = new AtomList();
		newAtomList.add(atom);
	    }
	}
	createAndAddResidue(creator, newAtomList, ResidueMode.CTER);
	sort();
    }


    public Chain(AtomList atomList, ResidueCreator creator, String chainName, Protein protein, MolecularSystem current,  MolecularSystem newms) {
	this(chainName, protein);
	AtomList newAtomList = null;
	boolean first = true;
	add(new Residue(new ResidueIdentifier(chainName,0))); //traditionally sequences start at 1
	if (atomList.size() == 0) throw new RuntimeException(" No Atoms in AtomList "+atomList.comment());
	for (Atom atom:atomList){
	    while (size()< atom.residueNumber()) {
		if (newAtomList == null) add(new Residue(new ResidueIdentifier(chainName,size())));
		else {
		    if (first) {
				createAndAddResidue(creator, newAtomList, ResidueMode.NTER);
				first = false;
		    }
		    else {
				createAndAddResidue(creator, newAtomList, ResidueMode.NORMAL);
			}
			newAtomList = null;
		}
	    }
	    if (size() == atom.residueNumber()) {
		if (newAtomList == null) newAtomList = new AtomList();
		newAtomList.add(atom);
	    }
	}
	createAndAddResidue(creator, newAtomList, ResidueMode.CTER);
	sort();
    }

    private static AtomList getCaAtomList(char c,Coordinates coordinates) {
	AtomType type = ResidueType.type(c).caType();
	if (type == AtomType.XXX) throw new RuntimeException("Problem in getCaAtomList"+
							     " cannot convert character \""+c+"\" to a residue type.");
	MolecularSystem saveMS = MolecularSystem.currentMolecularSystem();
	new MolecularSystem();
	Atom ca = new Atom("CA",null,type,coordinates,new Double(0));
	MolecularSystem.setCurrentMolecularSystem(saveMS);
	AtomList out = new AtomList();
	out.add(ca);
	return out;
    }
    private void createAndAddResidue(ResidueCreator creator, AtomList newAtomList, ResidueMode mode) {
	    ResidueIdentifier ri = new ResidueIdentifier(name(),size());
	    Residue residue = new Residue(ri);
	    try {
		residue = creator.create(newAtomList, ri, mode);
	    } catch (Exception ex) {
		System.out.println("Failed to create residue "+ri+"\n"+"from atom list:");
		newAtomList.print();
		System.out.println("dummy residue added");
	    }
	    add(residue);
	    residue.setResidueInAtoms(newAtomList);
	}

/*
	private void createAndAddResidue(ResidueCreator creator, int mode) {
		Residue residue = creator.create(newAtomList,size(),Residues.NORMAL);
	}
*/


    //------------------------------------------- methods -----------------------------------
    /**
     * Fetch a residue by its position in the list.
     **/
    public Residue residueAt(int index) {
	return get(index);
    }
    
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
	 * For users who need this signature.
 	 */
// 	public Residue residue(int residueNumber) {
// 		return residue(new Residue(residueNumber).ID());
// 	}

// Previous implementation:
/*
     public Residue residue(int residueNumber) {
         if (residueNumber >= size()) return null;
	 Residue out = get(residueNumber);
	 if (residueNumber < 0) 
	     throw new RuntimeException("Weird residueNumber "+residueNumber);
	 if (out.number != residueNumber) 
	     throw new RuntimeException("Resdiue "+out+" in position "+residueNumber+"\n"+
					"in "+this);
	 return out;
	 }
*/

// 	Residue key = new Residue(residueNumber);
// 	int index;
// 	try {
// 	    index = binarySearch(key);
// 	}
// 	catch (Exception ex) {
// 	    System.out.println("************ ERROR **************");
// 	    print(10);
// 	    ex.printStackTrace();
// 	    throw new RuntimeException("A problem with residueNumber "+residueNumber+"\n"+
// 				       "key = "+key+
// 				       "size = "+size()+"\n"+
// 				       "residues(residueNumber) = "+get(residueNumber)+"\n"+
// 				       ex+"\n");
// 	}
// 	if (index < 0) return null;
// 	return residueAt(index);

    /**
     * Extracts the residues accepted by the filter.
     **/
    public ResidueList filter(Filter residueFilter) {
        ResidueList newList = new ResidueList();
        for(Residue residue :this) {
            if (residueFilter.accept(residue))
            newList.add(residue);
        }
	return newList;
    }

    static class IsResidue implements Filter {
	public boolean accept(Object obj) {return (obj instanceof Residue);}
    } 



    
    public int numberOfNonDummyResidues() {
	int out = 0;
	for (Iterator residues = iterator(); residues.hasNext();) {
	    Residue residue = (Residue) residues.next();
	    if (!residue.dummy()) out++;
	}
	return out;
    }

	    


	public String toString() {
	String out = "Chain "+name+":\n";
	Iterator residues = iterator();
	Residue residue = (Residue) residues.next(); // remove the first dummy residue
	while (residues.hasNext()) {
	    residue = (Residue) residues.next();
	    if (residue.dummy()) out+=SequenceAlignmentCell.WILDCARD_CHAR;
	    else out+=residue.type().nameOneLetter();
	}
	return out;
    }
    
    public ResidueSequence sequence() {
	return new ResidueSequence(this, protein.name()+"_"+name());
    }
	
	public AtomList atoms() {return new AtomList(this);}

//     public void setChain(String chain) {
// //		assertSingleChain();
// 		Iterator residueIter = iterator();
// 		Residue residue;
// 		while (residueIter.hasNext()) {
// 			residue = (Residue) residueIter.next();
// 			residue.setChain(chain);
// 		}
// 	}

	/**
	 * Asserts that <code>this</code> chain is comprised of
	 * residues with the same chain string.
	 *
	 * @throws MeshiException if the assertion fails.
	 */
    public void assertChain() {
	for (Iterator residueIter = iterator();residueIter.hasNext();) {
	    Residue residue = (Residue) residueIter.next();
	    if (! residue.ID().chain().equals(name()))
		throw new MeshiException("Chain.assertChain -- " +
					 "It appears that this Chain "+this+" contains residue "+residue+"  of other chain.");
		}
    }

    public ResidueList missingResidues() {
	ResidueList out = new ResidueList();
	boolean nonDummyFound = false;
	for (Iterator residues = iterator(); residues.hasNext() & (! nonDummyFound);) {
	    Residue residue = (Residue) residues.next();
	    if (! residue.dummy()) nonDummyFound = true;
	}
	if (nonDummyFound) {
	    for (Iterator residues = iterator(); residues.hasNext() & (! nonDummyFound);) {
		Residue residue = (Residue) residues.next();
		if (residue.dummy()) out.add(residue);
	    }
	}
	return out;
    }
    public AtomList nowhereAtoms() {
	AtomList out = new AtomList();
	for (Iterator atoms = atoms().iterator();atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
	    if (atom.nowhere()) out.add(atom);
	}
	return out;
    }
    
    public Residue firstNonDummyResidue() {
	int index = 0;
	for (Iterator residues = iterator(); residues.hasNext();) {
	    Residue residue = (Residue) residues.next();
	    if (! residue.dummy()) {
		if (residue.number() != index) throw new RuntimeException("This is weird");
		return residue;
	    }
	    index++;
	}
	return null;
    }

}

