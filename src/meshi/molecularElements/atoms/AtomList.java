package meshi.molecularElements.atoms;
import meshi.molecularElements.*;
import meshi.util.*;
import meshi.util.file.*;
import meshi.util.filters.*;
import java.util.*;
import meshi.PDB.*;
import meshi.sequences.*; 
import meshi.parameters.*; 
/**
 * A list of Atoms
 **/
public class AtomList extends ArrayList<Atom> {
    private MeshiLineReader sourceFile = null;
    private MolecularSystem molecularSystem = null;
    public  MolecularSystem molecularSystem() {return molecularSystem;}
       //for SymmetryComplex
    public  void setMolecularSystem(MolecularSystem MolecularSystem) {this.molecularSystem = MolecularSystem;}
    private String comment = null;
   
    public AtomList() {
	super();
    }
    public AtomList(Atom atom) {
	super();
	add(atom);
    }
   public AtomList(int initialCapaciy) {
	super(initialCapaciy);
    }

    /**
     * An atom list based on a PDB file
     * In a multiple model file (such as NMR's) only model 1 is read.
     * Currently whenever there are two alternative atom positions present,
     * only position A is taken.
     **/
    public AtomList(MeshiLineReader file, Filter lineFilter){
	super();
	PdbLine line;
	boolean contRead = true;
	int maxUnknown = 1000;
	int nUnknown = 0;
	String[][] unknowns = new String[maxUnknown][2];
	sourceFile = file;
	if (!(file instanceof PdbReader)) 
	    throw new MeshiException("Cannot construct AtomList from file: "+file+"\n"+
			   "Currently only PDB formated files are used.");
	while (((line = ((PdbReader) file).readPdbLine()) != null) && contRead) {
	    if (line.isAModel() & (line.getModel() != 1))
		contRead = false;
	    if (lineFilter.accept(line)) {
		if (line.alternateLocation().equals("") || line.alternateLocation().trim().equals("A")) {
		    if (line.type().equals(AtomType.XXX)) {
			boolean found = false;
			for (int i = 0; (i < nUnknown) & (! found); i++) {
			    if (unknowns[i][0].equals(line.residueName()) & unknowns[i][1].equals(line.name())) 
				found = true;
			}
			if (! found) {
			    if (nUnknown >= maxUnknown) throw new RuntimeException("this is weird - too many unknown atom types");
			    System.out.println("Cannot identify Atom type for "+line);
			    unknowns[nUnknown][0] = line.residueName();
			    unknowns[nUnknown][1] = line.name();
			    nUnknown++;
			}
			    
		    } 
		    else {
			add(new Atom(line));
		    }
		}
		setComment(file.name());
	    }
	} 
	try {
	    file.close();
	}
	catch (Exception ex) { throw new RuntimeException("Cannot close "+file+"\n"+ex);}
    }
    /**
     * An atom list based on a PDB file
     **/
    public AtomList(PdbReader file) {
		this(file, new PdbLineATOM());
    }
    public AtomList(String dir, String fileName) {
		this(new PdbReader(dir, fileName));
    }
    public AtomList(String fileName, Filter lineFilter) {
		this(new PdbReader(fileName), lineFilter);
    }
    public AtomList(String fileName) {
		this(fileName, new PdbLineATOM());
    }
    /**
     * An atom list based on a residue list
     **/
    public  AtomList(ResidueList residues) {
	super();
	for (Residue residue:residues) {
	    if (residue.atoms() != null) {
		for (Atom atom:residue.atoms())
		    add(atom);
	    }
	}
    }
    /**
     * Extracts the list elements that pass the filter.
     **/
    public AtomList filter(Filter filter) {
	AtomList out = new AtomList();
	for (Atom atom:this)
	    if (filter.accept(atom)) out.add(atom);
	return out;
    }
	
    
    
    public boolean add(Atom atom){
	if (size() == 0) {
	    molecularSystem = atom.molecularSystem;
	}
	else {
	    if (atom.molecularSystem != molecularSystem)
		throw new RuntimeException("An attempt to add atom:\n"+
					   atom+"\n"+
					   "of "+atom.molecularSystem+" to an atom list of "+molecularSystem);
	}
	return(super.add(atom));
    }

    /**
     * RMS deviation between this list of atom and some other list.
     **/
    public Rms rms(AtomList otherList) {
	return new Rms(new AtomAlignment(this, otherList));
    }

    
    public double getRms(AtomList otherList) {
	return rms(otherList).getRms();
    }

    /**
     * Returns the atom in position <b>i</b> of the list.
     **/
    public Atom atomAt(int i) { return (Atom) get(i);}
    public static class IsAtom implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof Atom); 
        }
    }

    public int whereIs(Atom atom) {
	for (int i = 0; i < size(); i++) {
	    if (atomAt(i).equals(atom)) return i;
	}
	return -1;
    }

    /**
     * of gyration of the list atoms.
     * SQRT{ SIGMA_OVER_ALL_ATOMS([(x-X)^2+(y-Y)^2+(z-Z)^2]/N)}<br>
     * Where (x,y,z) are the atom coordinates and (X,Y,Z) are the 
     * coordinates of the center of mass.
     **/
    public double radius() {
              return  Math.sqrt(radius2() );
    }

    public double radius2() {
	double radius;
  	double dx, dy, dz, d2;
	double cmx, cmy, cmz; // center of mass x, y and z
	radius = 0;
	cmx = cmy = cmz = 0.0;
	int nWithCoordinates = 0;
	for (Atom atom:this) {
	  if (! atom.nowhere()) {
	    cmx += atom.x();
	    cmy += atom.y();
	    cmz += atom.z();
	    nWithCoordinates++;
	  }
	}
	cmx /= nWithCoordinates;
	cmy /= nWithCoordinates;
	cmz /= nWithCoordinates;
	for (Atom atom:this) {
	  if (! atom.nowhere()) {
		dx = cmx - atom.x();
		dy = cmy - atom.y();
		dz = cmz - atom.z();
		d2 = dx*dx + dy*dy + dz*dz;
		radius += d2;
	  }
	}
	return  radius / nWithCoordinates;
    }




    /**
     * Sets the parameter to be the residue of this atom.
     **/
    public void setResidue(Residue residue) {
	for (Atom atom:this)
	    atom.setResidue(residue);
    }

    public Atom[] toArrayOfAtoms() {
	return toArray(new Atom[size()]);
    }
    
    public MeshiLineReader sourceFile() { return sourceFile;}
    
    public boolean frozenAtomsExist() {
	for (Atom atom:this)
	    if (atom.frozen()) return true;
	return false;
    }

    public void freeze() {
	for (Atom atom:this)
	    if (atom.core.status == AtomStatus.NORMAL) atom.freeze();
    }
    public void freeze(Filter filter) {
	for (Atom atom:this) {
	    if (filter.accept(atom) & (atom.core.status == AtomStatus.NORMAL)) atom.freeze();
	}
    }

    public void defrost() {
	for (Atom atom:this)
	    if (atom.frozen()) atom.defrost();
    }
    public AtomList fullOccupancyFilter() {
	AtomList out = new AtomList();
	for (Atom atom:this)
	    if (atom.occupancy() > 0.99) out.add(atom);
	return out;
    }

    public AtomList noOXTFilter() {
	AtomList out = new AtomList();
	for (Atom atom:this)
	    if (atom.name().trim().compareTo("OXT") != 0) out.add(atom);
	return out;
    }

    public AtomList backbone() {
	return filter(new BackboneFilter());
    }

    public static class BackboneFilter implements Filter {
	public boolean accept(Object obj) {
	    Atom atom = (Atom) obj;
	    if (atom.name().equals("N")) return true; 
	    if (atom.name().equals("H")) return true; 
	    if (atom.name().equals("CA")) return true; 
	    if (atom.name().equals("CB")) return true; 
	    if (atom.name().equals("C")) return true;
	    if (atom.name().equals("O")) return true;
	    return false;
	}
    }

    public boolean sortable() {return false;}

    public String toString() { return "AtomList with "+size()+" atoms";}

    public AtomList CAFilter() {
	AtomList out = new AtomList();
 	out.sourceFile = sourceFile;
  	for (Atom atom:this)
           if (atom.name().trim().compareTo("CA") ==0) out.add(atom);
        return out;
    }
    
    public AtomList CAnoImageFilter() {
	AtomList out = new AtomList();
 	out.sourceFile = sourceFile;
  	for (Atom atom:this)
           if ((atom.name().trim().compareTo("CA") ==0) && (!atom.core.status().image())) out.add(atom);
        return out;
    }

	
     public AtomList frozenAtoms() {
       return filter(new IsFrozenFilter());
     }

     public AtomList defrostedAtoms() {
       return filter(new IsDefrostedFilter());
     }

     public static class IsDefrostedFilter implements Filter {
       public boolean accept(Object obj) {
           Atom atom = (Atom) obj;
           return !atom.frozen();
       }
     }

     public static class IsFrozenFilter implements Filter {
       public boolean accept(Object obj) {
           Atom atom = (Atom) obj;
           return atom.frozen();
       }
     }

     public static class NonHydrogen implements Filter {
       public boolean accept(Object obj) {
           Atom atom = (Atom) obj;
           return !(atom.isHydrogen());
       }
     }
    
    /**
     * Returns true if the parameter has the same length and the same composition in terms of atom names.
     * Current implementation, rather inefficient.
     **/
    public boolean sameAtoms(AtomList other) {
	if (size() != other.size()) return false;
	
	for (int iAtom = 0; iAtom < size(); iAtom++) {
	    boolean notFound =true;
	    Atom atom = get(iAtom);
	    for (int jAtom = 0; (jAtom < size()) & notFound; jAtom++) {
		Atom otherAtom = get(jAtom);
		if (atom.name().equals(otherAtom.name())) notFound = false;
	    }
	    if (notFound) return false;
	}
	return true;
    }
	    

	/**
	 * Asserts that <code>this</code> AtomList is comprised of
	 * atoms belonging to a single residue.
	 *
	 * @throws MeshiException if the assertion fails.
	 */
	public void assertSingleResidue() {
		if (! checkSingleResidue())
			throw new MeshiException("AtomList.assertSingleResidue -- " +
			"It appears that this AtomList contains more than one Residue.");
	}

	/**
	 * Checks if <code>this</code> AtomList is comprised of
	 * atoms belonging to a single residue.
	 */
	public boolean checkSingleResidue() {
	    if (size() == 0) return true;
		Iterator atomIter = iterator();
		Atom atom = (Atom) atomIter.next();
		Residue residueOfFirstAtom = atom.residue();

		while (atomIter.hasNext()) {
			atom = (Atom) atomIter.next();
			if (! Utils.equals(atom.residue(), residueOfFirstAtom))
				return false;
		}

		return true;
	}


	public AtomList[] splitToChains() {
		AtomList leftToSplit = this;
		ArrayList<AtomList> atomLists = new ArrayList<AtomList>();
		String currentChain;

		while (! leftToSplit.isEmpty()) {
			currentChain = leftToSplit.get(0).chain();

			atomLists.add(leftToSplit.filter(new AcceptChainFilter(currentChain)));
			leftToSplit = leftToSplit.filter(new RejectChainFilter(currentChain));
		}

		AtomList[] out = new AtomList[atomLists.size()];
		for (int i = 0; i < out.length; i++)
			out[i] = atomLists.get(i);

		return out;
	}

    public AtomList getChainAtoms(String chainLetter) { //A == 1
        return this.filter(new AcceptChainFilter(chainLetter));
    }

    public AtomList getSourceAtoms() {
        return this.splitToChains()[0];
    }


    public static class AcceptChainFilter implements Filter {
		public final String CHAIN;

		public AcceptChainFilter(String chain) {
			if (chain == null || chain.length() != 1)
				throw new MeshiException
					("ChainFilter must receive a one-letter String");

			this.CHAIN = chain;
		}

		public boolean accept(Object obj) {
			return (obj instanceof Atom) &&
					CHAIN.equals(((Atom) obj).chain());
		}
	}

	public static class RejectChainFilter implements Filter {
		public final String CHAIN;

		public RejectChainFilter(String chain) {
			if (chain == null || chain.length() != 1)
				throw new MeshiException
					("ChainFilter must receive a one-letter String");

			this.CHAIN = chain;
		}

		public boolean accept(Object obj) {
			return (obj instanceof Atom) &&
					! CHAIN.equals(((Atom) obj).chain());
		}
	}
    public AtomList somewhere() {
	AtomList out = new AtomList();
	for (Iterator atoms = iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
	    if (! atom.nowhere()) out.add(atom);
	}
	return out;
    }
    
    
    public void moveCMtoOrigin() {
  	double cmx = 0.0;
	double cmy = 0.0;
	double cmz = 0.0;
 	for (Atom atom:this) {
	    cmx += atom.x();
	    cmy += atom.y();
	    cmz += atom.z();
	}
	cmx /= size();
	cmy /= size();
	cmz /= size();
	for (Atom atom:this) {
	    atom.setXYZ(atom.x()-cmx,
			atom.y()-cmy,
			atom.z()-cmz);
	}
    }
    
    public boolean includesNowhereAtoms() {
	for (Iterator atoms = iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
	    if (atom.nowhere()) return true;
	}
	return false;
    }

    public void setComment(String s) {
	comment = s;
    }

    public void print() {
	for (Atom a:this) System.out.println(a);
    }
    public void print(MeshiWriter writer) {
	for (Atom a:this) writer.println(a);
    }

    public String comment() {return comment;}
}


