package meshi.molecularElements.ca;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;   
import meshi.util.*;
import meshi.util.string.*;
import meshi.util.filters.*;
import meshi.PDB.*;
import meshi.geometry.*;
import meshi.parameters.*;
import java.util.*;


public class CaResidue extends Residue{
    public static final CaResidue creator = new CaResidue("creator");
    public final Atom CA;    
    

   /**
     * <pre>
     * Use this constructor to instantiate a creator object..
     *
     **/
    public CaResidue(String name){
	super(name);
	CA = null;
    }

    public CaResidue(String name, ResidueType type, ResidueIdentifier id, AtomList atomList, ResidueMode mode) {
	super(name, type, id, getNewAtomList(atomList, type), mode);
	CA = atoms.get(0);
	CA.setResidue(this);
	head = tail = CA;
    }	
 
    public static AtomList getNewAtomList(AtomList oldAL, ResidueType type) {
	Atom old = find(oldAL,BBatom.CA);
	if (old == null) throw new RuntimeException("CA undefined");
	return new AtomList(new Atom("CA",null,type.caType(),new Coordinates(old),old.temperatureFactor()));
    }



    public int getContacts(ResidueList residueList, double contactDistance) {
	Iterator residues = residueList.iterator();
	Residue residue;
	double dis;
	int contacts = 0;
	while ((residue = (Residue) residues.next()) != null) {
	    if (!(residue.dummy())) {
		dis = (distance((CaResidue)residue)).distance();
		if (dis < contactDistance) contacts++;
	    }
	}
	return contacts;
    }

	
    public void setCoordinates(CaResidue from) {
	moveTo(from.CA.x(),
	       from.CA.y(),
	       from.CA.z());
    }
	
    public void setCoordinates(Atom from) {
	moveTo(from.x(),
	       from.y(),
	       from.z());
    }
	
    public  Distance distance(CaResidue other) {
	return new FreeDistance(CA,other.CA);
    }
    

    public void moveTo(double x, double y, double z) {
	CA.setXYZ(x,y,z);
    }

    public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode) {
	if ((atoms == null) || (atoms.size() == 0)) return new Residue(id);
	Atom atom = find(atoms,BBatom.CA);
	AtomType atomType = atom.type();
	ResidueType residueType = ResidueType.type(atomType);
	String residueName = residueType.nameThreeLetters();
        return new CaResidue(residueName, residueType,id, atoms, mode);
    } 
}
