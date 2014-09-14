package meshi.geometry;

import  meshi.molecularElements.atoms.*;
import java.util.*;
import meshi.util.*;

/**
 * A Distance Matrix specifically tailored for Calpha model. It is O(N2) but takes into account the geometrical constrains of CA model. If the chain
 * is not very long, it may run much faster than the assimptotically more efficient 
 *  
 **/
public class DistanceMatrixCaOnly extends DistanceMatrix {

    /*------------------------------------------ object variables --------------------------------------------*/
    private static int counter = 0;


    /*-------------------------------------------- constructors -----------------------------------------------*/
    public DistanceMatrixCaOnly(MolecularSystem molecularSystem, double rMax, int bondedListDepth) {
	super(molecularSystem);
    // Setting the contants used for calculating the reported distance.
    DistanceMatrix.rMax = rMax;
    this.bondedListDepth = bondedListDepth;
    reset();
    }
	public void refresh() {
	    counter = 0;
	}

    private void reset() {
	terminator.reset();
	rMax2 = rMax*rMax;
	int size = molecularSystem.size();
	
	// Initialize the distance matrix.
	matrix = new MatrixRow[size];
	for (int i = 0; i< size; i++) {
	    AtomCore atom = molecularSystem.get(i);
	    if (!atom.status().active()) {
		matrix[i] = null;
	    }
	    else matrix[i] = new MatrixRowCaOnly(atom, matrix, molecularSystem);
	}
	/**
	 * Every element of energyTermsDistanceLists
	 * can be specified in Creator of EnergyTerm needed it (for example HydrogenBondsCreator)
	 */
	energyTermsDistanceLists = new ArrayList<DistanceList>();
	
	// Generate bonded list.
	bondedList = getBondedList(molecularSystem, bondedListDepth, matrix);
	//bonded list.
 	try {
	    update();
	}
//	catch (UpdateableException e) {
    catch (Exception e) {
        System.out.println("Cannot update DistanceMatrixCaOnly");
        throw new RuntimeException("Cannot update DistanceMatrixCaOnly");
	}
    }
    
    
    /*------------------------------------------ update --------------------------------------------*/
    /**
     * Updates the distance matrix. 
     **/
    public void update(int numberOfUpdates) throws UpdateableException {
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    this.numberOfUpdates++;
	    update();
	}
	else if (numberOfUpdates != this.numberOfUpdates)
	    throw new RuntimeException("Something weird with DistanceMatrix.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }


    protected void update()  throws UpdateableException{
	if (terminator.dead()) throw new RuntimeException(terminator.message());

    if (nonBondedList == null) nonBondedList = new DistanceList(molecularSystem.size()*MatrixRow.CAPACITY);
    nonBondedList.clear();

    for (DistanceList dl:energyTermsDistanceLists)
	    dl.clear();

     //update bondedList
     for (Distance d:bondedList){
                 ((BondedDistance)d).update();
     }

    //build NBL    
    int length = matrix.length;
    for (int iRow = 0; iRow <length; iRow++) {
	if (matrix[iRow] != null)
	    ((MatrixRowCaOnly)matrix[iRow]).update(nonBondedList,energyTermsDistanceLists, counter);
    }
    counter++;

//        testNonBondedList();
    }


    /*------------------------------------------ other methods --------------------------------------------*/
    public void testNonBondedList(){
        Atom atom1, atom2;
        double dx,dy,dz,d2,dis;
         for (Distance distance:nonBondedList) {
             atom1 = distance.atom1();
             atom2 = distance.atom2();

        dx = atom1.x() - atom2.x();
		dy = atom1.y() - atom2.y();
		dz = atom1.z() - atom2.z();
		d2 = dx*dx + dy*dy + dz*dz;
             dis = Math.sqrt(d2);
             if (dis != distance.distance()){
                 System.out.println("Not Updated Distance! "+distance+" The real value is distance ="+dis);
            }

             DistanceList tmpRow = this.matrix[distance.atom1Number];
             boolean exist = false;
             for (Distance tmp:tmpRow){
             if (tmp == null) continue;
             if (tmp.atom2Number == distance.atom2Number){
                 exist = true;
                 if (tmp.distance() != distance.distance())
                    System.out.println("Mimatch between NBL and Matrix:"+distance+tmp);
                 break;
             }
            }
            if (!exist)
                    System.out.println("Mimatch between NBL and Matrix:"+distance+" does not exists in Matrix");

        }        
    }

    public String toString() {
    return ("DistanceMatrixCaOnly:\n"+
        "\t number of atoms \t"+molecularSystem.size()+
        "\t rMax \t"+rMax);
    }
}
