package meshi.symmetryComplex.molecularImageElements;

import meshi.util.*;
import java.util.*;
import meshi.molecularElements.*; // Unnecessary if in package molecularElements
import meshi.molecularElements.atoms.*;
import meshi.symmetryComplex.transformations.Transformation;
import static meshi.symmetryComplex.utils.GJUtils.getTransformations;
import static meshi.symmetryComplex.utils.GJUtils.RELATIVE_ORIENTATION;
import meshi.symmetryComplex.utils.GJFilters;
import meshi.sequences.*;
import meshi.parameters.SecondaryStructure;
import meshi.parameters.AtomType;
import meshi.geometry.Coordinates;
import meshi.geometry.Distance;
import meshi.geometry.FreeDistance;
import meshi.geometry.DistanceMatrix;

/**
 * Symmetric protein, each chain is a Protein object.
 * Chain A is "real", the others are images computed using Transformations.
 */

public class SymmetricComplex extends Protein implements Updateable {
    protected List<Transformation> transformations;
    public List<Transformation> transformations(){return transformations;}

    private int numberOfUpdates = 1;
    private String chainLetters;
    private static final String CHAIN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    private ResidueList imageResidues;
    public ResidueList imageResidues(){return imageResidues;}

    public SymmetricComplex(Transformation[] transformations, String chainLetters){
        super("SymmetricComplex");
        this.transformations = new ArrayList<Transformation>();

        for (Transformation transformation :transformations) {
            this.transformations.add(transformation);
        }
        this.chainLetters = chainLetters;
    }

    public SymmetricComplex(AtomList initialAtoms, ResidueCreator creator,
                            Transformation[] transformations){
         this(initialAtoms, creator, transformations, setConsecutiveChainLetters(transformations.length+1),false);
    }

    public SymmetricComplex(AtomList initialAtoms, ResidueCreator creator,
                            SequenceAlignment alignment,
                            Transformation[] transformations){
         this(initialAtoms, creator, transformations, setConsecutiveChainLetters(transformations.length+1),true);
         setSS(chains.get(0), alignment);
    }

    public SymmetricComplex(AtomList initialAtoms, ResidueCreator creator,
                            SequenceAlignment alignment,
                            Transformation[] transformations, String chainLetters){
         this(initialAtoms, creator, transformations, chainLetters,true);
         setSS(chains.get(0), alignment);
    }


    public SymmetricComplex(AtomList initialAtoms, ResidueCreator creator,
                            Transformation[] transformations, String chainLetters, boolean freeze){

        this(transformations, chainLetters);
//chain A from the sourceProtein
        new MolecularSystem();
        AtomList copyAtoms = new AtomList();//for new MolecularSystem()
        Atom copyAtom;
        for (Atom atom:initialAtoms){
                    copyAtom = new Atom(atom.name, atom.residue(), atom.type(), new Coordinates(atom), atom.temperatureFactor());
 //                   if (atom.frozen())   //does not work, as for Chain new Atom will be made without status
   //                             copyAtom.freeze();
                    copyAtoms.add(copyAtom);
        }

       new MolecularSystem();
       Chain chain = new Chain(copyAtoms, creator, "A", this);
       chains.add(chain);
       bonds = new AtomPairList(chains); //only real atoms yet
//add image atoms
       generateImageChains();

       residues = new ResidueList(chains);
       atoms = new AtomList(residues);
       ChainList imageChains =  chains.filter(new ImageChain.IsImageChain());// nonDummyResidues();
       imageResidues = new ResidueList(imageChains);
        if (freeze)
                atoms.freeze(GJFilters.frozenResiduesFilter);
    }
    
/*
    public SymmetricComplex(AtomList initialAtoms, ResidueCreator creator,
                            Transformation[] transformations, MolecularSystem current, MolecularSystem newms){
        this(transformations, setConsecutiveChainLetters(transformations.length+1));
//chain A from the sourceProtein
       Chain chain = new Chain(initialAtoms, creator, "A", this, current, newms);
       chains.add(chain);
       bonds = new AtomPairList(chains); //only real atoms yet
//add image atoms
       generateImageChains();

       residues = new ResidueList(chains);
       atoms = new AtomList(residues);

       ChainList imageChains =  chains.filter(new ImageChain.IsImageChain());// nonDummyResidues();
       imageResidues = new ResidueList(imageChains);

    }             */


     /*
    //toDo if I call it
 //  public void setTransformations(Transformation[] transformations, boolean [] chainsBooleans) {
//       setTransformations(transformations, setWorkingChainLetters(chainsBooleans));
//  }

   public void setTransformations(Transformation[] transformations, String chainLetters) {
        this.chainLetters = chainLetters;
        this.transformations = new ArrayList<Transformation>();
        for (Transformation transformation :transformations)
            this.transformations.add(transformation);

        if (chainLetters.length()-1 != transformations.length)
           throw new RuntimeException("Something weird in SymmetricComplex.setTransformation parameters");

        if (! chains.isEmpty())//todo
            generateImages();
    }
       */

    public void generateImageChains() {
        if ( chains.isEmpty() | transformations.isEmpty())
            throw new RuntimeException("Can't generate images: " +
                    "either there's no source or no transformations.");
         if (chainLetters.length()-1 != transformations.size())
            throw new RuntimeException("Try generate Images. Something weird with chains/transformations numbers in SymmetricComplex.");

        Chain source = getSource();          //chain 0
        ImageChain imageChain;

        for (int i = 1; i< chains.size();i++)
                        chains.get(i).clear();

        int i=1;
        for (Transformation transformation : transformations) {
            imageChain = new ImageChain(this,source,""+chainLetters.charAt(i),transformation);
//            imageChain = new ImageChain(this,source,""+chainLetters.charAt(i+1),transformation);
            chains.add(imageChain);
            i++;
        }
    }


    public Chain getSource() {
        return chains.get(0);
    }

    public List<Transformation> getTransformations() {
        return transformations;
    }

   //to update SymmetryComplex when all image have been built
    public void updateLocations() {

        for (Atom atom :atoms) {
            if (atom.core.status().image())
                    ((ImageAtom)atom).updateLocation();
        }
    }

    public void update(int numberOfUpdates) {
        if (numberOfUpdates == this.numberOfUpdates+1) {
            updateLocations();
            this.numberOfUpdates++;
        }
        else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with SymmetricComplex.update(int numberOfUpdates)\n"+
                    "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }

     public void resetNumberOfUpdates() {numberOfUpdates = 0;}

    //--------------------------------------------------------------------------------------------
    //method is the same to setSS in Protein, but this method set SS only to the part of residues from the first chain
    //
    public void setSS(ResidueList residues, SequenceAlignment alignment) {
        char prevSS = 'C';
        char ss;

    for (Iterator columns = alignment.iterator(); columns.hasNext();) {
        SequenceAlignmentColumn  column = (SequenceAlignmentColumn) columns.next();
        SequenceAlignmentCell cell = (SequenceAlignmentCell) column.cell(0);
        cell.addAttribute(column);
    }
    SequenceList tempList = new SequenceList(alignment);
    Sequence sequenceFromSS = tempList.get(0);
    SequenceAlignment newAlignment = SequenceAlignment.identityAlignment(chains.get(0).sequence(),sequenceFromSS);
    for(Residue residue:residues){
        if (!residue.dummy()) {
        SequenceAlignmentColumn alignmentColumn = (SequenceAlignmentColumn) newAlignment.getColumn(0,residue.ID().number());
               if (alignmentColumn.cell(1).gap()) ss = prevSS;
        else {
           SequenceAlignmentCell tempCell = (SequenceAlignmentCell) alignmentColumn.cell(1);
           SequenceAlignmentColumn tempColumn = (SequenceAlignmentColumn) tempCell.getAttribute();
           ss = tempColumn.getChar(1);
           prevSS = ss;
        }
           residue.setSecondaryStructure(SecondaryStructure.secondaryStructure(ss));
        }
    }
    }
//-----------------------------------------------------------------------------------------------

    public String chainLetters() {
        String letters = "";
        for (Chain chain : chains)
                letters += chain.name();
        return letters;
    }


    public static String setConsecutiveChainLetters(int size){
        return CHAIN_LETTERS.substring(0, size);
    }

    public void printChainsNames(){
        System.out.print("Chains in use: ");
        for (Chain chain : chains)
                System.out.print(chain.name()+" ");
        System.out.println();
    }

/*
//------------------------------------------------------------------------------------------------------------------------
  protected void cleanComplex( ) {
      if (chains != null) chains.clear();
      if (residues != null) residues.clear();
      if (atoms != null) atoms.clear();
      if (bonds != null) bonds.clear();
      if (angles != null) angles.clear();
      if (torsions != null) torsions.clear();
  }

  protected void cleanImageComplex() {
      Chain chain = this.getSource();
      chains.clear();
      chains.add(chain);

      if (residues != null) residues.clear();
      if (atoms != null) atoms.clear();

      if (angles != null) angles.clear();
      if (torsions != null) torsions.clear();
}
*/
}