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
 *           SG 
 **/
public class Cys extends ResidueExtendedAtoms {
    public final Atom SG;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"          SG \n";

    public Cys(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
   	super(ResidueType.CYS, atomList, id, mode);
	Object[] temp = new Object[1];
	SG = getAtom("SG",AtomType.CSG,atomList,this);
	atoms.add(SG);
	bonds.add(CB.bond(SG));
   }
    public String comment() {
	return COMMENT;
    }
}
