package meshi.symmetryComplex.geometryImage;

import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.AtomStatus;
import meshi.geometry.*;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Jan 7, 2009
 * Time: 2:41:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class MatrixRowNoImageAtoms extends MatrixRow {
    private static DistanceMode mode;
    private static boolean[] toNonBonded = new boolean[CAPACITY];

    public MatrixRowNoImageAtoms(AtomCore atom, MatrixRow[] matrix) {
    super(atom, matrix);
    }

    public void update(DistanceList nonBondedList, ArrayList<DistanceList> energyTermsDistanceLists) {
    Double x1, y1, z1, d2;
    x1 = atom.x();
    y1 = atom.y();
    z1 = atom.z();
    if (toNonBonded.length < size)
        toNonBonded = new boolean[2*capacity];
    for (int i = 0; i < size; i++) {
        toNonBonded[i] = false;
        distance = internalArray[i];
        if (distance != null) {
        mode = distance.mode();
        if (mode.mirror || mode.free)
            throw new RuntimeException("Distance mirrors and free distances have nothing to do here\n"+distance);
        if (!mode.frozen) {
            atom2 = distance.atom2;
            distance.dx = x1-atom2.x();
            distance.dy = y1-atom2.y();
            distance.dz = z1-atom2.z();
            d2 = distance.dx*distance.dx+distance.dy*distance.dy+distance.dz*distance.dz;
            if ((d2 < rMax2) || (mode.bonded)) {
            distance.distance = Math.sqrt(d2);
            distance.invDistance = 1/distance.distance;
            if ( (!mode.bonded) &&
                 (!( distance.atom1.status().image() && distance.atom2.status().image()))
                ){
                toNonBonded[i] = true;
                //nonBondedList.add(distance);
                if (mode == DistanceMode.INFINITE) {
                 mode = DistanceMode.NEW;
                distance.setMode(DistanceMode.NEW);
                for (DistanceList dl:energyTermsDistanceLists)
                    dl.add(distance);
                }
                else {
                distance.setMode(DistanceMode.NORMAL);
                }

            }
            }
            else if (d2 < DistanceMatrix.rMaxPlusBuffer2()) { //toDo to check value
            mode = DistanceMode.INFINITE;
            distance.setMode(DistanceMode.INFINITE);
            distance.distance = Distance.INFINITE_DISTANCE;
            distance.dx = distance.dy = distance.dz = 0;
            distance.invDistance = 1;
            }
            else  {
            internalArray[i] = null;
            distance.setMode(DistanceMode.DEAD); // just in case comebody is looking at this object
            }

        }
        }
    }

    for (int i = 0; i < size; i++)
        if (toNonBonded[i]) nonBondedList.add(internalArray[i]);
    }

    public void addCell(GridCell cell,DistanceList nonBondedList, ArrayList<DistanceList> energyTermsDistanceLists){
	double x = atom.x();
	double y = atom.y();
	double z = atom.z();
	double dx, dy, dz, d2, dis;
	Distance distance;
	double rMaxPlusBuffer2 = DistanceMatrix.rMaxPlusBuffer2();
	boolean found;
	for (AtomCore cellAtom:cell) {
	    int cellAtomNumber = cellAtom.number;
	    if (number > cellAtomNumber) {
		dx = x-cellAtom.x();
		dy = y-cellAtom.y();
		dz = z-cellAtom.z();
		d2 = dx*dx+dy*dy+dz*dz;
		if (d2 < rMaxPlusBuffer2) {
		    lastEmpty = -1;
		    found = false;
		    for (int i = 0; i < size; i++) {
			distance = internalArray[i];
			if (distance == null){
			    lastEmpty = i;
			}
			else {
			    atom2number = distance.atom2Number();
			    if (cellAtomNumber == atom2number) {
				found = true;
				    break;
			    }
			}
		    }
		    if (!found) {// needs to be inserted

			if (d2 < rMax2) dis = Math.sqrt(d2);
			else dis = Distance.INFINITE_DISTANCE;

			if ((atom.status() == AtomStatus.FROZEN)  & (cellAtom.status() == AtomStatus.FROZEN))
			    distance = new FrozenDistance(atom,cellAtom,dx,dy,dz,dis);
			else distance = new Distance(atom,cellAtom,dx,dy,dz,dis);
			if ((distance.distance < Distance.INFINITE_DISTANCE)&&
                (!( distance.atom1.status().image() && distance.atom2.status().image()))
                ){
			    distance.setMode(DistanceMode.NEW);
                nonBondedList.add(distance);
			    for (DistanceList dl:energyTermsDistanceLists)
				dl.add(distance);
			}
			else {
			    distance.setMode(DistanceMode.INFINITE);
			}
			if (lastEmpty == -1) add(distance);
			else internalArray[lastEmpty] = distance;
		    }
		}
	    }
	}
    }

}
