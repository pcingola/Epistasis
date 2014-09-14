package meshi.util;
import meshi.util.file.*;
import meshi.util.string.*;
import meshi.util.overlap.*;
import meshi.util.filters.*;
import meshi.sequences.*;
import meshi.geometry.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.ca.*;
import meshi.PDB.*;
import java.util.*;
import java.lang.*;
import java.io.*;


/**
 * A (hopefully) comfortable handle to the complexity of the Kabsch structural overlap algorithm[1,2]. 
 * The algorithm itself is implemented in meshi.overlap.Overlap.
 * 1. "A solution for the best rotation to relate two sets of vectors". By Wolfgang Kabsch, (1976) Acta Cryst. A 32: 922-923<BR>
 * 2. "A discussion of the solution for the best rotation to relate two sets<BR>
 *    of vectors". By Wolfgang Kabsch, Acta Cryst. (1978). A34, 827-828.
 * This class was rather drastically modified by Chen in version 1.45.
 **/
public class Rms implements KeyWords{

    private static boolean debug = true;
    private boolean alive = true;
    //data members
    private double[][] rotateMatrix;
    private double rms;
    private Coordinates centerOfMass0 = new Coordinates();
    private Coordinates centerOfMass1 = new Coordinates();
        

    //constructors

    public Rms(ResidueAlignment residueAlignment) {
	this(new AtomAlignment( residueAlignment));
    }
	
    public Rms(AtomAlignment atomAlignment) {
	if (atomAlignment.hasGaps()) throw new RuntimeException("Cannot calculate RMS for AtomAlignment with gaps");
	int size = atomAlignment.size();
	double[][] coor0 = new double[3][size];
	double[][] coor1 = new double[3][size];
	String comment0 = atomAlignment.get(0).cell(0).comment;
	String comment1 = atomAlignment.get(0).cell(1).comment;
	double centerOfMassX1 = 0, centerOfMassY1 = 0, centerOfMassZ1 = 0;
	double centerOfMassX0 = 0, centerOfMassY0 = 0, centerOfMassZ0 = 0;
	try {
	    for (int i = 0; i < size; i++){
		Atom atom0 = atomAlignment.atomAt(i,0);
		Atom atom1 = atomAlignment.atomAt(i,1);
		coor0[0][i] = atom0.x();
		coor0[1][i] = atom0.y();
		coor0[2][i] = atom0.z();
		coor1[0][i] = atom1.x();
		coor1[1][i] = atom1.y();
		coor1[2][i] = atom1.z();
		
		centerOfMassX0 += atom0.x();
		centerOfMassY0 += atom0.y();
		centerOfMassZ0 += atom0.z();
		centerOfMassX1 += atom1.x();
		centerOfMassY1 += atom1.y();
		centerOfMassZ1 += atom1.z();
	    }
	    centerOfMassX0 /= size;
	    centerOfMassY0 /= size;
	    centerOfMassZ0 /= size;
	    centerOfMassX1 /= size;
	    centerOfMassY1 /= size;
	    centerOfMassZ1 /= size;

	    centerOfMass0.setXYZ(centerOfMassX0,centerOfMassY0,centerOfMassZ0);
	    centerOfMass1.setXYZ(centerOfMassX1,centerOfMassY1,centerOfMassZ1);
	    
	    Overlap overlap = new Overlap(coor0, coor1, size, comment0, comment1);
	    rms = overlap.rms();
	    rotateMatrix = overlap.rotationMatrix();
	}
	catch (Exception e) {
	    throw new MeshiException("Rms Error:\n"+
				  "comparing "+comment0+"\n"+
				  "with\n"+comment1+"\n"+
				  "reproted problem"+
				  e);
	}
    }
	

    private Coordinates centerOfMass0() { return centerOfMass0;}
    private Coordinates centerOfMass1() { return centerOfMass1;}

    public static double superimpose(Protein protein0, Protein protein1, ResidueAlignment residueAlignment){
	AtomAlignment atomAlignment = new AtomAlignment(residueAlignment);
	 Rms rms = new Rms(atomAlignment);
	 Coordinates centerOfMass0 = rms.centerOfMass0();
	 Coordinates centerOfMass1 = rms.centerOfMass1();
	 for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
	     Atom atom = (Atom) atoms.next();
	     atom.addToX(0-centerOfMass1.x());
	     atom.addToY(0-centerOfMass1.y());
	     atom.addToZ(0-centerOfMass1.z());
	 }
	 double[][] rotateMatrix = rms.getMatrix();
	 for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
	     Atom atom = (Atom) atoms.next();
	     double x = atom.x();
	     double y = atom.y();
	     double z = atom.z();
	     atom.setXYZ(rotateMatrix[0][0]*x + rotateMatrix[0][1]*y + rotateMatrix[0][2]*z,
			 rotateMatrix[1][0]*x + rotateMatrix[1][1]*y + rotateMatrix[1][2]*z,
			 rotateMatrix[2][0]*x + rotateMatrix[2][1]*y + rotateMatrix[2][2]*z);
	 }
	 for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
	     Atom atom = (Atom) atoms.next();
	     atom.addToX(centerOfMass0.x());
	     atom.addToY(centerOfMass0.y());
	     atom.addToZ(centerOfMass0.z());
	 }
	 return rms.getRms();
    }
	 	
	




    //methods

    public double[][] getMatrix(){
	if (alive) return rotateMatrix;
	else return null;
    }

  
    public double getRms(){
    	if (alive) return rms;
	else return -1;
    }

    public String toString() {
	return (new Double(getRms())).toString();
    }

    public Coordinates getCenterOfMass0() {return centerOfMass0;}
    public Coordinates getCenterOfMass1() {return centerOfMass1;}

    //------------------------------------------------------------------ utility methods --------------------------------------------------------------

     public  static final double[] NO_GDT = {-1, -1, -1, -1, -1};

    public static double rms(Protein protein0, CommandList commands) {
	return rms(protein0, commands, new KolDichfin()); // KolDichfin - a filter that accept all
    }

    public static double rms(Protein protein0, CommandList commands, Filter filter) {
	AtomAlignment atomAlignment = getAtomAlignment(protein0, commands, filter);
	if (atomAlignment == null) return -1;
	Rms rms = new Rms(atomAlignment);
	return rms.getRms();
    }

    public static double rms(Protein protein0, Protein protein1) {
	ResidueAlignment alignment = new ResidueAlignment(protein0.chain(), protein0.name(), protein1.chain(), protein1.name());
	return rms(alignment);
    }

    public static double[] gdt(Protein protein0, Protein protein1) {
	return gdt( new ResidueAlignment(protein0.chain(), protein0.name(),protein1.chain(),protein1.name()), protein0.atoms().CAFilter().size());
    }
    
    public static double[] gdt (ResidueAlignment residueAlignment, int refLength) {
	return gdt(residueAlignment, new KolDichfin(), refLength);
    }

    public static double rms (ResidueAlignment residueAlignment) {
	AtomAlignment atomAlignment = new AtomAlignment(residueAlignment, new KolDichfin());
	if (atomAlignment == null) return -1;
	Rms rms = new Rms(atomAlignment);
	return rms.getRms();
    }

    public static double[] gdt(Protein protein0, CommandList commands) {
	return gdt(protein0, commands, new KolDichfin()); // KolDichfin - a filter that accept all
    }

    public static double[] gdt(Protein model, CommandList commands, Filter filter) {
	//	try {
		CommandList alignmentCommands = commands.firstWordFilter(SUPERIMPOSE);
		String referenceFileName = alignmentCommands.secondWord(REFERENCE).thirdWord();
		if (referenceFileName.equals(NONE.key)) {
			 	System.out.println("No reference structure");
				return NO_GDT;
		}
		Protein reference = Utils.getProtein(commands, SUPERIMPOSE, REFERENCE, CaResidue.creator,Utils.defaultExceptionHandler);
		int refLength = reference.atoms().CAFilter().size();
		ResidueAlignment residueAlignment = new ResidueAlignment(reference.chain(),reference.name(),model.chain(),model.name());
		return gdt(residueAlignment, filter, refLength);
		//	}
		//tch (Exception ex) {System.out.println("\n\n"+"Failed to calculate gdt\n"+ex);}
		//turn NO_GDT;
    }

    public static double[] gdt (ResidueAlignment residueAlignment, Filter filter, int refLength) {
	//y {
	        AtomAlignment atomAlignment = new AtomAlignment(residueAlignment,filter);
		if (atomAlignment == null ) return NO_GDT;
		AtomList referenceAtoms = atomAlignment.atomList(0);
		AtomList modelAtoms     = atomAlignment.atomList(1);

	        AtomList newReferenceAtoms = Utils.duplicateInAnewMolecularSystem(referenceAtoms);
	        AtomList newModelAtoms     = Utils.duplicateInAnewMolecularSystem(modelAtoms);

	        return GDTcalculator.gdt(newReferenceAtoms, newModelAtoms, refLength);
		//
		//tch (Exception ex) {System.out.println("\n\n"+"Failed to calculate gdt\n"+ex);}
		//turn NO_GDT;
    }
	    


    private static AtomAlignment getAtomAlignment(Protein protein0, CommandList commands, Filter filter) {
	CommandList alignmentCommands = commands.firstWordFilter(SUPERIMPOSE);
	String referenceFileName = alignmentCommands.secondWord(REFERENCE).thirdWord();
	if (referenceFileName.equals(NONE.key)) {
	    System.out.println("No reference structure");
	    return null;
	}
	else {
	    Protein protein1 = Utils.getProtein(commands, SUPERIMPOSE, REFERENCE, CaResidue.creator, Utils.defaultExceptionHandler);
	    if (debug) System.out.println("protein1 = "+protein1);
	    String mode = alignmentCommands.secondWord(MODE).thirdWord();
	    if (!mode.equals(ALL_CA.key)) 
		throw new RuntimeException("Currently the only implemented mode of "+SUPERIMPOSE+"\n"+
					   "is "+ALL_CA);
	    System.out.println("Generating ResidueAlignment of \n"+protein1.sequence()+"\n"+protein0.sequence());
	    ResidueAlignment residueAlignment = new ResidueAlignment(protein1.chain(), protein1.name(),protein0.chain(),protein0.name());
	    Utils.print(residueAlignment);
	    return  new AtomAlignment(residueAlignment, filter);
	}
    }
}
