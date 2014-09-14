package meshi.molecularElements.extendedAtoms;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.geometry.putH.*;
import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.angle.*;
/**
 *<pre>
 *ARGININ
 *
 *              O
 *              |
 *     N - CA - C...n
 *     |   |
 *     H   CB                 
 *         |                  
 *         CG - CD - NE - CZ - NH1
 *                        |
 *                        NH2
  **/
public class Arg extends ResidueExtendedAtoms {
    public final Atom CG, CD, NE, HE, CZ, NH1, NH2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"              O\n"+
	"              |\n"+
	"     N - CA - C...n\n"+
	"     |   |\n"+
	"     H   CB\n"+                 
	"         |\n"+                  
	"         CG - CD - NE - CZ - NH1\n"+
	"                   |    |\n"+
	"                   HE   NH2\n";


  public Arg(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.ARG, atomList, id, mode);
	
	atoms.add(CG = getAtom("CG",AtomType.RCG, atomList, this));
	atoms.add(CD = getAtom("CD",AtomType.RCD, atomList, this));
	atoms.add(NE = getAtom("NE",AtomType.RNE, atomList, this));
	atoms.add(HE = getAtom("HE",AtomType.RHE, atomList, this));
	atoms.add(CZ = getAtom("CZ",AtomType.RCZ, atomList, this));
	atoms.add(NH1 = getAtom("NH1",AtomType.RNH, atomList, this));
	atoms.add(NH2 = getAtom("NH2",AtomType.RNH, atomList, this));
	bonds.add(CB.bond(CG));
	bonds.add(CG.bond(CD));
	bonds.add(CD.bond(NE));
	bonds.add(NE.bond(CZ));
	bonds.add(NE.bond(HE));
	bonds.add(CZ.bond(NH1));
	bonds.add(CZ.bond(NH2));
    }

    public void addHydrogens(BondParametersList bondParametersList, AngleParametersList angleParametersList) {
	try {
	    super.addHydrogens(bondParametersList, angleParametersList); 
	    if (HE.nowhere() && 
		(!CD.nowhere()) &&
		(!NE.nowhere())&&
		(!CZ.nowhere())) {
		PutHpos.pos(HE,bondParametersList, angleParametersList);
	    }
	    //	    if ((NE != null) && HE.nowhere()) PutHpos.pos(HE,bondParametersList, angleParametersList);
	} catch (NotEnoughBoundAtomsException ex) {HE.resetCoordinates();}
    }

    public String comment() {
	return COMMENT;
    }
}
