package meshi.symmetryComplex.molecularImageElements;

import java.util.*;
import meshi.molecularElements.*; // Unnecessary if in package molecularElements
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomStatus;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.util.MeshiException;

/**
 * A list of ImageAtoms
 **/
public class ImageAtomList extends AtomList  {

    protected Transformation transformation;

    public ImageAtomList(AtomList sourceAtomList, Transformation transformation) {
        super();
        this.transformation = transformation;
        if (sourceAtomList != null) {
            Iterator atoms = sourceAtomList.iterator();
            Atom atom;
            while (atoms.hasNext()) {
                atom = (Atom) atoms.next();
                add(new ImageAtom(atom, transformation));  
            }

        }
    }

    public ImageAtomList(Residue residue, Transformation transformation) {
        this(residue.atoms(), transformation);
    }

    public  ImageAtomList(ImageChain residues, Transformation transformation) {
        super();
        this.transformation = transformation;
        Iterator resIter = residues.iterator();
        ImageResidue residue;
        Iterator atomsIter;
        ImageAtom atom;
        while (resIter.hasNext()) {
            residue = (ImageResidue) resIter.next();
            atomsIter = residue.atoms().iterator();
            while (atomsIter.hasNext()) {
               atom = (ImageAtom) atomsIter.next();
               if (atom.core.status() != AtomStatus.IMAGE)
                   throw new MeshiException ("This method should add image atoms only.");
                    atom.freeze();
               add(atom);
            }
        }
    }
}