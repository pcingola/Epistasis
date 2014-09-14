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
 *       H   CB
 *           |
 *           CG - OD2
 *           | (-)
 *           OD1
 **/
public class Asp extends ResidueExtendedAtoms {
    public final Atom CG, OD1, OD2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"           CG - OD2\n"+
	"           | \n"+
	"           OD1\n";

    public Asp(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.ASP, atomList, id, mode);

	atoms.add(CG = getAtom("CG",AtomType.DCG, atomList, this));
	atoms.add(OD1 = getAtom("OD1",AtomType.DOD, atomList, this));
	atoms.add(OD2 = getAtom("OD2",AtomType.DOD, atomList, this));
	bonds.add(CB.bond(CG));
	bonds.add(CG.bond(OD2));
	bonds.add(CG.bond(OD1));
    }
    public String comment() {
	return COMMENT;
    }
}

