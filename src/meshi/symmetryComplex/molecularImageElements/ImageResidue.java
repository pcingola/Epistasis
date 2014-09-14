package meshi.symmetryComplex.molecularImageElements;

import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.Atom;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.parameters.ResidueType;
import meshi.parameters.ResidueMode;

public class ImageResidue extends Residue {

    private Residue source;

    public ImageResidue(Residue source, String chainName, Transformation transformation) {
        super (source.name, source.type, new ResidueIdentifier(chainName,source.number()),
                new ImageAtomList(source.atoms(), transformation),
                source.mode);
        this.source = source;
//        this.dummy = source.dummy();//todo
        setSecondaryStructure(source.secondaryStructure());
    }

    public Residue source(){return source;}
}
