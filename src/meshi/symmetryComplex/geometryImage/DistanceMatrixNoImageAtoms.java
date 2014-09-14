package meshi.symmetryComplex.geometryImage;

import java.util.ArrayList;

import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.util.UpdateableException;
import meshi.geometry.*;

public class DistanceMatrixNoImageAtoms extends DistanceMatrix {
    private Indicator  indicatorToUpdateHB;

     public DistanceMatrixNoImageAtoms(MolecularSystem molecularSystem, double rMax,
                                       double buffer, double edge, int bondedListDepth, ResidueList imageResidues) {
    super(molecularSystem);
	// Setting the contants used for calculating the reported distance.
	nonBondedList =  new DistanceList(molecularSystem.size()*10);
	this.buffer = buffer;
	DistanceMatrix.rMax = rMax;
	this.edge = edge;
	this.bondedListDepth = bondedListDepth;
	reset();

// add C-O, H-N bond distances for Ami hbAngle energyTerm
      Atom cAtom, oAtom, hAtom, nAtom;
      Distance distance;
      double dx, dy, dz, d;
         if (imageResidues == null) return;
             for (Residue residue: imageResidues){
                 oAtom = residue.carbonylO();    //in Meshi oAtom.number > cAtom.number always
                 cAtom = residue.carbonylC();
                 if ((cAtom != null) && ( oAtom != null)) {
                     dx = oAtom.x() - cAtom.x();
                     dy = oAtom.y() - cAtom.y();
                     dz = oAtom.z() - cAtom.z();
                      d = Math.sqrt(dx*dx+dy*dy+dz*dz);
                     distance = new BondedDistance(oAtom, cAtom, dx, dy, dz, d);

                      if ((matrix[oAtom.number()] == null) || ( matrix[cAtom.number()] == null))
                        throw new RuntimeException("Too small  matrix for the imageResidueList");
                      if (oAtom.number() > cAtom.number())  //always true
                           matrix[oAtom.number()].add(distance);
                        else
                          throw new RuntimeException("Weird order of atoms in residue (oAtom before cAtom) "+ residue);
                 }

                 nAtom = residue.amideN();
                 hAtom = residue.amideH();  //in Meshi nAtom.number > hAtom.number always
                 if ((hAtom != null) && ( nAtom != null)) {
                     dx = nAtom.x() - hAtom.x();
                     dy = nAtom.y() - hAtom.y();
                     dz = nAtom.z() - hAtom.z();
                      d = Math.sqrt(dx*dx+dy*dy+dz*dz);
                     distance = new BondedDistance(nAtom, hAtom, dx, dy, dz, d);

                      if ((matrix[nAtom.number()] == null) || ( matrix[hAtom.number()] == null))
                        throw new RuntimeException("Too small  matrix for this imageResidueList");
                      if (nAtom.number() > hAtom.number())  //always
                           matrix[nAtom.number()].add(distance);
                        else
                          throw new RuntimeException("Weird order of atoms in residue (nAtom before hAtom) "+ residue);
                 }
             }
    }

    private void reset() {
    terminator.reset();
    rMax2 = rMax*rMax;
    rMaxPlusBuffer  = rMax+buffer;
    rMaxPlusBuffer2  = rMaxPlusBuffer*rMaxPlusBuffer;
    bufferOneThirdSqr = buffer*buffer/9;
    int size = molecularSystem.size();
    indicatorToUpdateHB = new Indicator();

    // Initialize the distance matrix.
    matrix = new MatrixRowNoImageAtoms[size];
    for (int i = 0; i< size; i++) {
        AtomCore atom = molecularSystem.get(i);
        if (!(atom.status().activeOrImage())){   
        matrix[i] = null;
        }
        else matrix[i] = new MatrixRowNoImageAtoms(atom, matrix);
    }
    /**
     * Every element of energyTermsDistanceLists
     * can be specified in Creator of EnergyTerm needed it (for example HydrogenBondsCreator)
     */
    energyTermsDistanceLists = new ArrayList<DistanceList>();

    // Generate bonded list.
    bondedList = getBondedList(molecularSystem, bondedListDepth, matrix);
    //bonded list.

    // We adjust the grid cell edge size so that the memory requirements of the
    // grid will not result in memory failure.
    try {
        while ((grid = new Grid(molecularSystem, edge, DEFAULT_EDGE(rMax,buffer))).failToBuild()) {
        if (edge >= 100000) throw new RuntimeException("Very weird");
           edge*=2;
        }


        update();
    }
    catch (UpdateableException e) {
        throw new RuntimeException(" Cannot create grid. Apparently the atoms a spread over"+
                       " a very large volume");
    }
    }

    protected void update()  throws UpdateableException{
	if (nonBondedList == null) nonBondedList = new DistanceList(molecularSystem.size()* MatrixRow.CAPACITY);
	nonBondedList.clear();

	for (DistanceList dl:energyTermsDistanceLists)
	    dl.clear();
	if (terminator.dead()) throw new RuntimeException(terminator.message());
	int length = matrix.length;
	for (int iRow = 0; iRow <length; iRow++)
	    if (matrix[iRow] != null) {
		matrix[iRow].update(nonBondedList,energyTermsDistanceLists);
	    }
	int size = molecularSystem.size();

        if (! grid.build()) throw new UpdateableException();
        for (int iatom = 0; iatom < size; iatom++) {
	    MatrixRow row = matrix[iatom];
	    if (row != null) {
		GridCell gridCell = grid.getCell(row.atom);
		row.addCell(gridCell,nonBondedList,energyTermsDistanceLists);
	    }
	}
                                                /*
	    for (Distance distance :nonBondedList){
            if (distance.atom1.status().image())
                System.out.println("ImageAtom in NBL in Distance:"+distance+distance.atom1());
            if (distance.atom2.status().image())
                System.out.println("ImageAtom in NBL in Distance:"+distance+distance.atom2());
            if (distance.atom1.status().image() && distance.atom2.status().image())
                throw new RuntimeException("Both atoms in NBL distance are image!"+distance);

        }                                         */
        //  testNonBondedList();
    }
    
}
