package meshi.symmetryComplex.proteinsGJ;

import java.io.IOException;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ResidueCreator;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.ca.CaResidue;
import meshi.molecularElements.extendedAtoms.BackboneResidue;
import meshi.util.CommandList;
import meshi.util.Command;
import meshi.util.Utils;
import meshi.symmetryComplex.utils.GJFilters;
import static meshi.symmetryComplex.utils.GJUtils.getTransformations;
import static meshi.symmetryComplex.utils.GJUtils.RELATIVE_ORIENTATION;
import static meshi.symmetryComplex.utils.GJUtils.whichChainsInZone;
import static meshi.symmetryComplex.utils.GJUtils.getPartialTransformations;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.geometry.fragments.DresserFullFrag;
import meshi.PDB.PdbLineATOM;


public class BackboneCompleterAndMinimizer extends GJModeller  {

    public BackboneCompleterAndMinimizer(String sequence,
                                         CommandList commands,
                                         String outputFileMask) throws IOException {
        super(sequence, commands, outputFileMask);
    }

    public BackboneCompleterAndMinimizer(String tmCasPdbLocation, String sequence,
                                         CommandList commands,
                                         String outputFileMask, Boolean generateAllData) throws IOException {
//CA
        super(tmCasPdbLocation, sequence, commands,
                outputFileMask);
        if (!toContinue("CA")) return;
        caRelax();
        //update();
       if (!toContinue("CA")) {
           System.out.println("Crosses or twists have been determined after caRelaxation.");
           return;
       }
        restoreAll(creator);
        generateFiles(gj, creator);
        timeStamp();
        generateLoopData("CA");
        
//BACKBONE
        backboneMinimization();

        if (!toContinue("CB")) return;
        generateFiles(gj, creator);
        timeStamp();

        if (generateAllData){
            restoreAll(creator);
            generateData("CB");
            generateLoopData("CB");
            //generateCysMatrixFile(gj.atoms(), "CB");
        }
    }

 protected void backboneMinimization() throws IOException{
        backboneMinimization(null);
}

    protected void backboneMinimization(String templatefileName) throws IOException{
    // complete backbone atoms
    creator = new BackboneResidue("creator");
    new MolecularSystem();
    Protein bbp;
    if (templatefileName ==  null) {
        bbp = new Protein(gj.getSource().atoms(),creator);

        new MolecularSystem();
        Command command = commands.firstWord(PARAMETERS_DIRECTORY);
        String FragmentsFileName  = command.secondWord();
        command = commands.firstWord(DRESSER_FRAGMENTS);
        FragmentsFileName  = FragmentsFileName+"/"+command.secondWord();
        DresserFullFrag dresser = new DresserFullFrag(FragmentsFileName);
        AtomList tempAtomList = dresser.dressProt(gj.getSource().atoms());
        Utils.assignBackboneCoordinates(bbp.atoms(),tempAtomList);
    }
       else     {
                bbp = new Protein(templatefileName,new PdbLineATOM(),creator);
    }


      // freeze TM Cas
      //positionTMCas(tmCasProtein, bbp, predictedResiduesFilter);
    //bbp.atoms().freeze(GJFilters.frozenResiduesFilter);  //now it is in the SymmetryComplex

    new MolecularSystem();
    SymmetricComplex gjAll = new SymmetricComplex(bbp.atoms().getSourceAtoms(),creator, alignment, getTransformations(RELATIVE_ORIENTATION,numberOfChains));

    generateFiles(gjAll,creator);
    timeStamp();

    gj = makeActualSymmetryComplex(gjAll, creator);
    // minimize backbone model
     if (templatefileName ==  null)
         backboneRelax();
   }

}
