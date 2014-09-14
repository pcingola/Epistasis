package meshi.symmetryComplex.molecularImageElements;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.util.filters.Filter;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Dec 23, 2008
 * Time: 11:44:13 AM
 * To change this template use File | Settings | File Templates.
 */
public class ImageChain extends Chain {
    private Transformation transformation;
    private Chain source;
    public Chain source(){return source;}
    public Transformation transformation() {return transformation;}

    public ImageChain(Protein protein, Chain source, String chainName, Transformation transformation){
        super(chainName,protein);     //todo check that this chain is empty
        this.source = source;
        this.transformation = transformation;
        for (Residue residue :source) {
            add(new ImageResidue(residue, chainName, transformation));
        }
    }

   public static class IsImageChain implements Filter {
     public boolean accept(Object obj) {return (obj instanceof ImageChain);}
     }

}
