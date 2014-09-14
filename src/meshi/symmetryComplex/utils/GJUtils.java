package meshi.symmetryComplex.utils;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Dec 23, 2008
 * Time: 2:56:18 PM
 * To change this template use File | Settings | File Templates.
 */

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.geometry.Coordinates;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.symmetryComplex.transformations.SingleAxisRotationTransformation;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;
import meshi.util.Utils;

/**
 * This class contains utils for work on gap junctions. So far all of them are
 * static -- use static import.
 *
 */
public class GJUtils {
    public static final double RELATIVE_ORIENTATION = 43.0;
    /**
     * Residue number ranges of the gap junction's parts.
     */
    //public static final int[]
//	                      //  N_TER = { 1, 18 },
                      //      M1 = { 19, 40 },
                       //     E1 = { 41, 75 },
                        //  M2 = { 76, 96 },
                         //   IN = { 97, 130 };
                          //  M3 = { 131, 152 },
                       //     E2 = { 153, 188 },
                           // M4 = { 189, 209 },
                  //          C_TER = { 210, 999 };
/*/
    public static final int[]
	                        N_TER = { 1, 18 },
	                        M1 = { 19, 39 },
	                        E1 = { 40, 75 },
	                        M2 = { 76, 96 },
	                        IN = { 97, 132 },
	                        M3 = { 133, 152 },
	                        E2 = { 153, 194 },
	                        M4 = { 195, 211 },
	                        C_TER = { 210, 999 };
//*/
    //public static int[][] E = { E1, E2 };
    //public static int[][] tmResiduesRanges = { M1, M2, M3, M4 };

    public static HydrogenBondsCreator hBondsCreator = new HydrogenBondsCreator();

    /**
     * Coordinates of transmembranal C-Alphas from Sar-el's work were deposited
     * in the PDB under the entry 1txh. This function takes a protein built from
     * those Cas and a whole monomer, and positions and freezes the
     * transmembranal Cas according to 1txh.
     *
     */
    public static void positionAndFreezeCas(final Protein partialCas, final Protein monomer) {
        Iterator monomerCasIterator;
        final Iterator partialCasIterator = partialCas.atoms().iterator();
        Atom monomerAtom, partialCasAtom;
        int partialCasNumber;
        outer: while (partialCasIterator.hasNext()) {
            partialCasAtom = (Atom) partialCasIterator.next();
            partialCasNumber = partialCasAtom.residueNumber();
            monomerCasIterator = monomer.atoms().CAFilter().iterator();
            while (monomerCasIterator.hasNext()) {
                monomerAtom = (Atom) monomerCasIterator.next();
                if (partialCasNumber == monomerAtom.residueNumber()) {
                    monomerAtom.setXYZ(partialCasAtom.x(), partialCasAtom.y(),
                            partialCasAtom.z());
                    monomerAtom.freeze();
                    continue outer;
                }
            }
        }
    }



    public static Transformation[] getTransformations(final double relativeOrientation, int numberOfChains) {
        int halfNumberOfChains = numberOfChains/2;

        final Transformation[] ans = new Transformation[numberOfChains - 1];
        final Transformation apposingTransformation = new SingleAxisRotationTransformation(
                "x", 180.0);
        Transformation firstHexamerTrans;
        Transformation secondHexamerTrans;

        ans[halfNumberOfChains-1] = apposingTransformation
                .compose(new SingleAxisRotationTransformation("z",
                        relativeOrientation));
        for (int i = 1; i <= halfNumberOfChains-1; i++) {
            firstHexamerTrans = new SingleAxisRotationTransformation("z", -(360/halfNumberOfChains)*i);
            secondHexamerTrans = new SingleAxisRotationTransformation("z", (360/halfNumberOfChains)*i + relativeOrientation);
            ans[i - 1] = firstHexamerTrans;
            ans[i + (halfNumberOfChains-1)] = apposingTransformation.compose(secondHexamerTrans);
        }
        return ans;
    }
/*
	public static Transformation[] getHalfTransformations(double relativeOrientation) {
		Transformation[] all = getTransformations(relativeOrientation);
		Transformation[] ans = new Transformation[BoundariesMap.halfNumberOfChains];
		ans[0] = all[0];
		ans[1] = all[4];
		ans[2] = all[5];
		ans[3] = all[6];
		ans[4] = all[9];
		ans[5] = all[10];
		return ans;
	}
*/
    /**
     * Coordinates of transmembranal C-Alphas from Sar-el's work were deposited
     * in the PDB under the entry 1txh. This function takes a protein built from
     * those Cas and a whole monomer, and positions the transmembranal Cas
     * according to 1txh.
     */
    public static void positionTMCas(final Protein tmCasProtein, final Protein monomer ) {
        final Iterator tmCasIterator = tmCasProtein.atoms().iterator();
        Atom monomerCa, tmCa;
        while (tmCasIterator.hasNext()) {
            tmCa = (Atom) tmCasIterator.next();
            monomerCa = monomer.getAtom(tmCa.residueName(), tmCa.residue().ID(), "CA");//5.29 tmCa.residueNumber();
            monomerCa.setXYZ(tmCa.x(), tmCa.y(), tmCa.z());
        }
    }

    /**
     * Adds atoms to the protein, builds a complex using the transformations and
     * writes the atom list to a file with the given name. Useful for
     * Ca-proteins, which SPDBV refuses to open.
     */                                         /*
    public static void spdbv_ify(final AtomList atoms, final Transformation[] transformations,
                                 final String filename) throws IOException {
        spdbv_ify(atoms, transformations).print(new MeshiWriter(filename));
    }
                                                               //*/
    /**
     * Adds atoms to the protein, builds a complex using the transformations and
     * returns the atom list. Useful for Ca-proteins, which SPDBV refuses to
     * open.
     */                                                  /*

    public static AtomList spdbv_ify(final AtomList atoms,
                                     final Transformation[] transformations) {
        final Protein p = new Protein(new AtomList(atoms), (new ResidueExtendedAtomsCreator()));   //5.29

        int maxNumberOfRandomCoordinatesPerResidue = 100;
        for (Iterator residues = p.residues().iterator();residues.hasNext();){
            Residue residue = (Residue) residues.next();
            if (!residue.dummy()) {
            boolean OK = Utils.assignRandomCoordinates(residue, maxNumberOfRandomCoordinatesPerResidue);
                    if (! OK )
                        throw new RuntimeException("Can not assign RandomCoordinates");
                }
        }

        final SymmetricComplex s = new SymmetricComplex(p, new ResidueExtendedAtomsCreator(), null,transformations);
        return s.atoms();
    }
                                                                    //*/
    /**
     * Calculates the average distance between atoms from loops E1 and E2. If
     * the given atom list contains more than one chain, only atoms from the
     * first chain will be taken into account.
     */
    /*public static double getAverageDistanceBetweenLoops(AtomList atoms) {
     atoms = atoms.splitToChains()[0];
     final Atom[] e1Array = atoms.filter(new ResidueRangeFilter(E1))
             .toArrayOfAtoms(), e2Array = atoms.filter(
             new ResidueRangeFilter(E2)).toArrayOfAtoms();
     double sum = 0.0;
     for (final Atom e1Atom : e1Array)
         for (final Atom e2Atom : e2Array)
             sum += e1Atom.distanceFrom(e2Atom);
     return sum / (e1Array.length * e2Array.length);
 }
    */


    public static boolean[] whichChainsInZone(double zone, AtomList atoms) {
        AtomList[] chains = atoms.splitToChains();
        boolean[] ans = new boolean[chains.length-1];
        Atom[] chainAArray = chains[0].toArrayOfAtoms(), chainArray;
        Atom[][] noChainAArrays = new Atom[chains.length - 1][];
        for (int i = 1; i < chains.length; i++)
            noChainAArrays[i - 1] = chains[i].toArrayOfAtoms();

    chainSearch:
        for (int chain = 0; chain < noChainAArrays.length; chain++) {
            ans[chain] = false;
            chainArray = noChainAArrays[chain];
            for (int i = 0; i < chainAArray.length; i++) {
                for (int j = 0; j < chainArray.length; j++)
                    if (quickDistanceFrom(chainAArray[i].core.x(),chainAArray[i].core.y(),chainAArray[i].core.z(),         //5.29
                                        chainArray[j].core.x(),chainArray[j].core.y(),chainArray[j].core.z()) <= zone) {
                        ans[chain] = true;
                        continue chainSearch;
                    }
            }
        }
        return ans;
    }

    public static String whichChainsInZone(ChainList chainList, double zone) {
       // AtomList[] chains = atoms.splitToChains();
        String actualChainsNames = "A";
        Atom[] chainAArray = chainList.get(0).atoms().toArrayOfAtoms(),
                     chainArray;
        Atom[][] noChainAArrays = new Atom[chainList.size() - 1][];
        for (int i = 1; i < chainList.size(); i++)
            noChainAArrays[i - 1] = chainList.get(i).atoms().toArrayOfAtoms();

    chainSearch:
        for (int chain = 0; chain < noChainAArrays.length; chain++) {
            chainArray = noChainAArrays[chain];
            for (int i = 0; i < chainAArray.length; i++) {
                for (int j = 0; j < chainArray.length; j++)
                    if (quickDistanceFrom(chainAArray[i].core.x(),chainAArray[i].core.y(),chainAArray[i].core.z(),         //5.29
                                        chainArray[j].core.x(),chainArray[j].core.y(),chainArray[j].core.z()) <= zone) {
                        actualChainsNames += (chainList.get(chain+1)).name();
                        continue chainSearch;
                    }
            }
        }

        return actualChainsNames;
    }


    public static Transformation[] getPartialTransformations(double relativeOrientation, int numberOfChains, boolean[] whichTrans) {
        Transformation[] all = getTransformations(relativeOrientation, numberOfChains);
        if (whichTrans.length != all.length)
            throw new RuntimeException("Booleans array length ("+whichTrans.length+
                    ") is different than allTrans array length ("+all.length+").");
        List<Transformation> ans = new Vector<Transformation>();
        for (int i = 0; i < all.length; i++) {
            if (whichTrans[i])
                ans.add(all[i]);
        }
        return ans.toArray(new Transformation[] {});
    }

    public static double quickDistanceFrom(double x1, double y1, double z1, double x2, double y2, double z2){//Coordinates fromCoors, Coordinates toCoors) {
        double dx = x2 - x1,
        dy = y2 - y1,
        dz = z2 - z1;
        return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }

        public static double quickDistanceFrom(Coordinates fromCoors, Coordinates toCoors) {
        double dx = toCoors.x() - fromCoors.x(),
        dy = toCoors.y() - fromCoors.y(),
        dz = toCoors.z() - fromCoors.z();
        return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }
}
