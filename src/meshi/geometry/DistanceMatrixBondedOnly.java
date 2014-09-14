package meshi.geometry;
import  meshi.molecularElements.*;
import  meshi.molecularElements.atoms.*;
import  meshi.parameters.*;
import java.util.*;
import meshi.util.filters.*;
import meshi.util.*;
/**
 
 **/
public class DistanceMatrixBondedOnly extends DistanceMatrix {
	
	
    private static double dx,dy,dz,d2;
    private int numberOfUpdates = 0;
    private Distance[] distanceArray;

    /*-------------------------------------------- constructors -----------------------------------------------*/
    public DistanceMatrixBondedOnly(MolecularSystem molecularSystem, int bondedListDepth) {
	super(molecularSystem);
	this.bondedListDepth = bondedListDepth;
	reset();
    }	
	
    private void reset() {
	bondedList = getBondedList(molecularSystem, bondedListDepth);
	bondedList.sort();
	distanceArray = new Distance[bondedList.size()];
	for (int i = 0; i < bondedList.size(); i++)
	    distanceArray[i] = (Distance) bondedList.get(i);
	Arrays.sort(distanceArray, new DistanceComparator());
    }

    public static DistanceList getBondedList(MolecularSystem molecularSystem, int depth) {
	AtomList bonded;
	Atom atom, bondedAtom;
	Iterator bondedAtoms;
	DistanceList out;
	int size = molecularSystem.size(); 
	AtomPairList tempList = new AtomPairList();
	for (int iatom = 0; iatom < size; iatom++) {
	    atom = molecularSystem.get(iatom).atom;
	    bonded = getBonded(atom, depth);
	    bondedAtoms = bonded.iterator();
	    while ((bondedAtom = (Atom) bondedAtoms.next()) != null) 
		tempList.add(new AtomPair(atom, bondedAtom));
	}
	tempList.sort();
	// This instantiate the Distance objects associated with bonded atom pairs. 
	// These pairs will be ignored during the non-bonded list update.	
	AtomPair atomPair;	
	Distance distance;
	out = new DistanceList(100);
	Iterator atomPairs = tempList.iterator();
	while ((atomPair = (AtomPair) atomPairs.next()) != null) {
		dx = atomPair.atom1().x() - atomPair.atom2().x();
		dy = atomPair.atom1().y() - atomPair.atom2().y();
		dz = atomPair.atom1().z() - atomPair.atom2().z();
		d2 = dx*dx+dy*dy+dz*dz;
		if (atomPair.atom1Number() > atomPair.atom2Number())
		    distance = new BondedDistance(atomPair.atom1(),atomPair.atom2(),dx,dy,dz,Math.sqrt(d2));
	    else 
		distance = new BondedDistance(atomPair.atom2(),atomPair.atom1(),-dx,-dy,-dz,Math.sqrt(d2));
 	    out.add(distance);
	    out.add(new DistanceMirror(distance));
	}
	return out;
    }

     public static AtomList getBonded(Atom atom, int depth) {
	AtomList out = new AtomList();
	getBonded(atom, depth, out,atom.number());
	return out;
    }
    public static void  getBonded(Atom atom, int depth, AtomList out, int rootNumber) {        
	if (depth == 0) return;
	Iterator atoms = atom.bonded().iterator();	
	Atom bondedAtom;
	while ((bondedAtom = (Atom) atoms.next()) != null) {
	    if ((rootNumber < bondedAtom.number()) &
		(! out.contains(bondedAtom)))
		out.add(bondedAtom);
	    getBonded(bondedAtom, depth-1, out, rootNumber);
	}
	
    }

    /**
     * Updates the distance matrix. 
     **/
    public void update(int numberOfUpdates) throws UpdateableException {
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    this.numberOfUpdates++;
            //if (numberOfUpdates % 2 == 1)
       //    if ((numberOfUpdates < 50) ||  (numberOfUpdates % 20 == 1))
	                          update();
            
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with DistanceMatrix.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }				       

    protected void update()  throws UpdateableException{
	for (Distance distance:distanceArray) 
	    if (!distance.mode.mirror) {
		distance.dx = distance.atom1().x() - distance.atom2().x(); 
		distance.dy = distance.atom1().y() - distance.atom2().y(); 
		distance.dz = distance.atom1().z() - distance.atom2().z(); 
		d2 = distance.dx*distance.dx+distance.dy*distance.dy+distance.dz*distance.dz;
		distance.distance = Math.sqrt(d2);
	    }
    }
   

    /** 
     * Returns the Distance object of the parameters.
     **/
    public Distance distance(Atom atom1, Atom atom2) {
	return distance(atom1.number(), atom2.number());
    }
   public Distance distance(AtomPair atomPair) { 
	return distance(atomPair.largeNumber(), atomPair.smallNumber());
    }
    public Distance distance(int atom1Number, int atom2Number) {
	Distance distance;
	int high = distanceArray.length-1;
	int low = 0;
	int middle;
	while  (low <= high )   {
	    middle = (low+high)>>1;
	    distance =  distanceArray[middle];
	    //System.out.println("yyyyy "+atom1Number+" "+atom2Number+" | "+low+" "+high+" "+middle+" | "+distance);
	    if (distance.atom1Number > atom1Number) {high = middle - 1;}
	    else {
		if (distance.atom1Number == atom1Number) {
		    if (distance.atom2Number > atom2Number) {high = middle - 1;}
		    else {
			if (distance.atom2Number == atom2Number) return distance;
			else {low = middle + 1;}
		    }
		}
		else {low = middle + 1;}
	    }
	}
	throw new RuntimeException("distance "+atom1Number+" to "+atom2Number+" not found");
    }

    private class DistanceComparator implements Comparator {
	public int compare (Object obj1, Object obj2) {
	    Distance distance1 = (Distance) obj1;
	    Distance distance2 = (Distance) obj2;
	    if (distance1.atom1Number < distance2.atom1Number) return -1;
	    if (distance1.atom1Number > distance2.atom1Number) return 1;
	    if (distance1.atom2Number < distance2.atom2Number) return -1;
	    if (distance1.atom2Number > distance2.atom2Number) return 1;
	    return 0;
	}
    }




    public String toString() {
	return ("DistanceMatrixBondedOnly:\n"+
		"\t number of atoms \t"+molecularSystem.size()+
		"\t rMax \t"+rMax+
		"\t buffer\t"+buffer); 
    }


    public void testNonBondedList() {
	System.out.println("No non-bonded-list to test");
    }
}
