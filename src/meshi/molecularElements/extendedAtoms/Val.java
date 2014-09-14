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
 *       H   CB - CG1
 *           |
 *           CG2
 *
 **/
public class Val extends ResidueExtendedAtoms {
    public final Atom CG1, CG2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
 "                O\n"+
 "                |\n"+
 "       N - CA - C...n\n"+
 "       |   |\n"+
 "       H   CB - CG1\n"+
 "           |\n"+
 "           CG2\n";


    public Val(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.VAL, atomList, id, mode);
	atoms.add(CG1 = getAtom("CG1",AtomType.VCG1, atomList, this));
	atoms.add(CG2 = getAtom("CG2",AtomType.VCG2, atomList, this));
	bonds.add(CB.bond(CG1));
	bonds.add(CB.bond(CG2));
    }
    public String comment() {
	return COMMENT;
    }
}
