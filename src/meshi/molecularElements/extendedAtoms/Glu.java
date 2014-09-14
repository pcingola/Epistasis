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
 *           CG - CD - OE2
 *                | (-)
 *                OE1
 **/
public class Glu extends ResidueExtendedAtoms {
    public final Atom CG, CD, OE1, OE2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"           CG - CD - OE2\n"+
	"                |\n"+
	"                OE1\n";

    public Glu(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.GLU, atomList, id, mode);
	atoms.add(CG = getAtom("CG",AtomType.ECG, atomList, this));
	atoms.add(CD = getAtom("CD",AtomType.ECD, atomList, this));
	atoms.add(OE1 = getAtom("OE1",AtomType.EOE, atomList, this));
	atoms.add(OE2 = getAtom("OE2",AtomType.EOE, atomList, this));
	bonds.add(CB.bond(CG));
	bonds.add(CG.bond(CD));
	bonds.add(CD.bond(OE1));
	bonds.add(CD.bond(OE2));
    }
    public String comment() {
	return COMMENT;
    }
}
