package meshi.molecularElements.extendedAtoms;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.parameters.*;
/**
 *<pre>
 *
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB - CG1 - CD1
 *           |
 *           CG2
 *
 **/
public class Ile extends ResidueExtendedAtoms {
    public final Atom CG1, CD1, CG2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB - CG1 - CD1\n"+
	"           |\n"+
	"           CG2\n";
 
   public Ile(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.ILE, atomList, id, mode);
 
	atoms.add(CG1 = getAtom("CG1",AtomType.ICG1, atomList, this));
	atoms.add(CD1 = getAtom("CD1",AtomType.ICD, atomList, this));
	atoms.add(CG2 = getAtom("CG2",AtomType.ICG2, atomList, this));
	bonds.add(CB.bond(CG1));
	bonds.add(CG1.bond(CD1));
	bonds.add(CB.bond(CG2));
    }
    public String comment() {
	return COMMENT;
    }
}
