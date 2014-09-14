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
 *       CD  CB
 *       \  /
 *        CG
 **/
public class Pro extends ResidueExtendedAtoms {
    public final Atom CG, CD;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       CD  CB\n"+
	"       \\  /\n"+
	"        CG\n";


  public Pro(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.PRO, atomList, id, mode);
	atoms.add(CG = getAtom("CG",AtomType.PCG, atomList, this));
	atoms.add(CD = getAtom("CD",AtomType.PCD, atomList, this));
	bonds.add(CB.bond(CG));
	bonds.add(CG.bond(CD));
	bonds.add(CD.bond(N));
    }
    public String comment() {
	return COMMENT;
    }
}
