package meshi.molecularElements.extendedAtoms;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.parameters.*;
/**
 *<pre>
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB             
 *           |              
 *           CG - CD - CE - NZ
 *
 **/
public class Lys extends ResidueExtendedAtoms {
    public final Atom CG, CD, CE, NZ;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+             
	"           |\n"+              
	"           CG - CD - CE - NZ\n";

  public Lys(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.LYS, atomList, id, mode);
 
	atoms.add(CG = getAtom("CG",AtomType.KCG, atomList, this));
	atoms.add(CD = getAtom("CD",AtomType.KCD, atomList, this));
	atoms.add(CE = getAtom("CE",AtomType.KCE, atomList, this));
	atoms.add(NZ = getAtom("NZ",AtomType.KNZ, atomList, this));
	bonds.add(CB.bond(CG));
	bonds.add(CG.bond(CD));
	bonds.add(CD.bond(CE));
	bonds.add(CE.bond(NZ));
    }
    public String comment() {
	return COMMENT;
    }
}
