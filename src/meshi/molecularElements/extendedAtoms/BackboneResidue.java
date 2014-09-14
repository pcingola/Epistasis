package meshi.molecularElements.extendedAtoms;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.geometry.*;
/**
 *<pre>
 * Alanin from Levitt, JMB 168:592 (1983) table 2.
 *           CB  O
 *           |   |
 *      N - CA - C...n
 *      |
 *      H
 **/
public class BackboneResidue extends ResidueExtendedAtoms implements ResidueCreator {
    public static final String COMMENT = "General backbone residue.";
    public static final BackboneResidue creator = new BackboneResidue("creator");

    public BackboneResidue(String name) {
   	super(name);
    }

	public BackboneResidue(ResidueType type, AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
   	super(type, atomList, id, mode);
    }

    public String comment() {
	return COMMENT;
    }
    public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode) {
	if ((atoms == null) || (atoms.size() == 0)) return new Residue(id);
	Atom atom = find(atoms,BBatom.CA);
	AtomType atomType = atom.type();
	ResidueType residueType = ResidueType.type(atomType);
        return new BackboneResidue(residueType, atoms, id, mode);
    } 
}
