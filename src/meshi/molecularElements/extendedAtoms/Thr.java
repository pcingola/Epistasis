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
 *           |
 *           CB
 *          /  \
 *        OG1  CG2
 **/
public class Thr extends ResidueExtendedAtoms{
    public final Atom OG1, CG2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"           |\n"+
	"           CB\n"+
	"          /  \\\n"+
	"        OG1  CG2\n";

 

  public Thr(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.THR, atomList, id, mode);

	atoms.add(CG2 = getAtom("CG2",AtomType.TCG, atomList, this));
	atoms.add(OG1 = getAtom("OG1",AtomType.TOG, atomList, this));
	bonds.add(CB.bond(OG1));
 	bonds.add(CB.bond(CG2));
    }
    public String comment() {
	return COMMENT;
    }
}
