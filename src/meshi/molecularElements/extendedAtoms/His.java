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
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB
 *           |
 *           CG - CD2 - NE2 - HE2
 *           |           |
 *     HD1 - ND1 ------ CE1
 **/
public class His extends  ResidueExtendedAtoms{
    public final Atom CG, CD2, NE2, HE2, CE1, ND1, HD1;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"               O\n"+
	"               |\n"+
	"      N - CA - C...n\n"+
	"      |   |\n"+
	"      H   CB\n"+
	"          |\n"+
	"          CG - CD2 - NE2 - HE2\n"+
	"          |           |\n"+
	"    HD1 - ND1 ------ CE1\n";

 
    public His(AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(ResidueType.HIS, atomList, id, mode);

	atoms.add(CG = getAtom("CG",AtomType.HCG, atomList, this));
	atoms.add(ND1 = getAtom("ND1",AtomType.HND, atomList, this));
	atoms.add(HD1 = getAtom("HD1",AtomType.HHD, atomList, this));
	atoms.add(CE1 = getAtom("CE1",AtomType.HCE, atomList, this));
	atoms.add(NE2 = getAtom("NE2",AtomType.HNE, atomList, this));
	atoms.add(HE2 = getAtom("HE2",AtomType.HHE, atomList, this));
	atoms.add(CD2 = getAtom("CD2",AtomType.HCD, atomList, this));
	bonds.add(CB.bond(CG));
 	bonds.add(CG.bond(ND1));
 	bonds.add(ND1.bond(HD1));
 	bonds.add(ND1.bond(CE1));
 	bonds.add(CG.bond(CD2));
 	bonds.add(CD2.bond(NE2));
 	bonds.add(CE1.bond(NE2));
 	bonds.add(NE2.bond(HE2));
   }
    public String comment() {
	return COMMENT;
    }
    public void addHydrogens(BondParametersList bondParametersList, AngleParametersList angleParametersList) {
	super.addHydrogens(bondParametersList, angleParametersList);
	try {
	    if (HD1.nowhere() && 
		(!CG.nowhere()) &&
		(!ND1.nowhere())&&
		(!CE1.nowhere())) {
		PutHpos.pos(HD1,bondParametersList, angleParametersList);
	    }
	    if (HE2.nowhere() && 
		(!CD2.nowhere()) &&
		(!NE2.nowhere())&&
		(!CE1.nowhere())) {
		PutHpos.pos(HE2,bondParametersList, angleParametersList);
	    }
	} catch (NotEnoughBoundAtomsException ex) {HD1.resetCoordinates();HE2.resetCoordinates();}
    }
}
