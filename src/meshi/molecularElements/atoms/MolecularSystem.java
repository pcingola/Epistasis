package meshi.molecularElements.atoms;
import meshi.util.*;
import meshi.parameters.*;
import java.util.*;

public class MolecularSystem extends ArrayList<AtomCore> {
    public final int ID;
    private static int numberOfMolecularSystems = 0;
    private static MolecularSystem currentMolecularSystem = new MolecularSystem();
    private int numberOfAtoms = 0;

    public MolecularSystem() {
	super();
	ID = numberOfMolecularSystems;
	numberOfMolecularSystems++;
	currentMolecularSystem = this;
    }

    public static MolecularSystem currentMolecularSystem() { return currentMolecularSystem;}
    public static void setCurrentMolecularSystem(MolecularSystem ms){ 
	currentMolecularSystem = ms;
    }


    public String toString() {
	return "Molecular System "+ID+"     ";
    }
   
    protected AtomCore  createAtomCore(Atom atom, AtomType type, AtomStatus status, double  x, double y, double z) {
	AtomCore newAtomCore = new AtomCore(atom,          type,            status, size(), x,        y,        z);
	add(newAtomCore);
	return newAtomCore;
    }
}
