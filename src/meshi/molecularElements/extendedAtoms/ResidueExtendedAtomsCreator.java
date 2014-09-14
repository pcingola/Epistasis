package meshi.molecularElements.extendedAtoms;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.geometry.*;
import meshi.util.*;
import java.util.*;
public class ResidueExtendedAtomsCreator implements ResidueCreator {
    public static final ResidueExtendedAtomsCreator creator = new ResidueExtendedAtomsCreator();
    public ResidueExtendedAtomsCreator() {}
    
   public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode){
	Atom atom = atoms.atomAt(0);
	String name;
	ResidueType type;
	if (atom == null) throw new RuntimeException("This is weird "+id);
	if (atom.residue() != null) {
	    name = atom.residue().name;
	    type = ResidueType.type(name);
	}
        else {
	    type =  ResidueType.type(atom.type());
	}
	    
	ResidueExtendedAtoms out = null;
	if (type == ResidueType.ALA) out = new Ala(atoms, id, mode);
	if (type == ResidueType.CYS) out = new Cys(atoms, id, mode);
	if (type == ResidueType.ASP) out = new Asp(atoms, id, mode);
	if (type == ResidueType.GLU) out = new Glu(atoms, id, mode);
	if (type == ResidueType.PHE) out = new Phe(atoms, id, mode);
	if (type == ResidueType.GLY) out = new Gly(atoms, id, mode);
	if (type == ResidueType.HIS) out = new His(atoms, id, mode);
	if (type == ResidueType.ILE) out = new Ile(atoms, id, mode);
	if (type == ResidueType.LYS) out = new Lys(atoms, id, mode);
	if (type == ResidueType.LEU) out = new Leu(atoms, id, mode);
	if (type == ResidueType.MET) out = new Met(atoms, id, mode);
	if (type == ResidueType.ASN) out = new Asn(atoms, id, mode);
	if (type == ResidueType.PRO) out = new Pro(atoms, id, mode);
	if (type == ResidueType.GLN) out = new Gln(atoms, id, mode);
	if (type == ResidueType.ARG) out = new Arg(atoms, id, mode);
	if (type == ResidueType.SER) out = new Ser(atoms, id, mode);
	if (type == ResidueType.THR) out = new Thr(atoms, id, mode);
	if (type == ResidueType.VAL) out = new Val(atoms, id, mode);
	if (type == ResidueType.TRP) out = new Trp(atoms, id, mode);
	if (type == ResidueType.TYR) out = new Tyr(atoms, id, mode);
	return out;
    }
}
