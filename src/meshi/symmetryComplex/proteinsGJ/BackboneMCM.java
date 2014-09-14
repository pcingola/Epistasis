package meshi.symmetryComplex.proteinsGJ;

import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.Command;
import meshi.molecularElements.ca.CaResidue;
import meshi.molecularElements.ResidueCreator;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.BackboneResidue;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.CooperativeAtomicPairwisePMFSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.inflate.simpleInflate.SimpleInflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.simpleInflate.SimpleInflate;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativePropensityCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.applications.prediction.homology.CaAtomFilter;
import meshi.symmetryComplex.geometryImage.DistanceMatrixNoImageAtoms;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.fragments.DresserFullFrag;
import meshi.optimizers.MCM;
import meshi.optimizers.Perturbation;
import meshi.optimizers.Minimizer;
import meshi.optimizers.TemperatureGenerator;
import meshi.symmetryComplex.utils.GJFilters;
import static meshi.symmetryComplex.utils.GJUtils.getTransformations;
import static meshi.symmetryComplex.utils.GJUtils.RELATIVE_ORIENTATION;
import static meshi.symmetryComplex.utils.GJUtils.whichChainsInZone;
import static meshi.symmetryComplex.utils.GJUtils.getPartialTransformations;
import meshi.symmetryComplex.energy.SymmetryCreator;
import meshi.symmetryComplex.energy.cylinder.CylinderCreator;
import meshi.symmetryComplex.energy.edmEnergy.EDMEnergyCreator;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;
import meshi.symmetryComplex.transformations.Transformation;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 20/02/2008
 * Time: 14:48:03
 * To change this template use File | Settings | File Templates.
 */
public class BackboneMCM extends BackboneCompleterAndMinimizer  {

   
    public BackboneMCM  (String tmCasPdbLocation, String sequence, CommandList commands,
                         String outputFileMask) throws IOException {
          super(tmCasPdbLocation, sequence, commands,
                outputFileMask, false);

          backboneMCM();

          generateFiles(gj, creator);
          timeStamp();

          restoreAll(creator);
          generateData("CB-MCM");
          generateLoopData("CB-MCM");
          //generateCysMatrixFile(gj.atoms(), "CB-MCM");
        }

  protected void backboneMCM() throws IOException {
      System.out.println("MCM");
       gj.printChainsNames();

 //     SolvateCBCreator solvateCBCreator = new SolvateCBCreator();
      SimpleInflateCreator inflateCreator = new SimpleInflateCreator();
      HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
      AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
      RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
     CompositePropensityCreator propensityCreator =  new CompositePropensityCreator();

      EnergyCreator[] creators = {
              new SymmetryCreator(),
              new BondCreator(),
              new AngleCreator(),
              new PlaneCreator(),
              new OutOfPlaneCreator(),

//              solvateCBCreator,
              ramachCreator,
              new CooperativeRamachandranCreator(ramachCreator),
              propensityCreator,
              new CooperativePropensityCreator(propensityCreator),
//              new ExcludedVolCreator(),
              summaCreator,
              new CooperativeAtomicPairwisePMFSummaCreator(summaCreator),
//              new SolvateCreatorHBforMinimization(1.0 , 1.0, 1.0, 1.0,commands),
              hydrogenBondsCreator,
              new HydrogenBondsPairsCreator(hydrogenBondsCreator),
              new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
              new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
               //new TetherCreator(),

            new EDMEnergyCreator(GJFilters.chainA_Ca_loops_Atom_Filter),
   //         getCysCBetasDistConstraintsCreator(),
            new CylinderCreator(commands, "innerBarrel", new GJFilters.AndFilter(new CaAtomFilter(), innerBarrelFilter)),
            new CylinderCreator(commands, "outerBarrel", new GJFilters.AndFilter(new CaAtomFilter(), outerBarrelFilter)),
            inflateCreator,
          };
//          DistanceMatrixNoImageAtoms distanceMatrix =
//                  new DistanceMatrixNoImageAtoms(gj, solvateCBCreator.getParametersList(commands));
//      DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj.molecularSystem(), gj.fullMolecularSystem(),
  //                        6.5, 2.0, DistanceMatrix.DEFAULT_EDGE(6.5,2.0), 4, gj.actualImageResidues());
      DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj.atoms().molecularSystem(),
                              5.5, 2.0, DistanceMatrix.DEFAULT_EDGE(5.5,2.0), 4, gj.imageResidues());

          energy = new TotalEnergy(gj, distanceMatrix, creators, commands);

  //       energy.test();

          System.gc();
          try {
              MCM mcm = getMCM(gj, energy, commands, (SimpleInflate) inflateCreator.term());
              mcm.run(null);
          }
          catch (Exception ex) {
              generateFiles(gj, creator);
              energy.test();
              throw new RuntimeException(ex);}

      }

  private MCM  getMCM(Protein model, TotalEnergy energy, CommandList commands, SimpleInflate inflate) {
      System.out.println("xxxxxxxx "+inflate);
  Perturbation         perturbation         = new Perturbation(energy,commands, inflate, model,0);
  Minimizer            minimizer            = Utils.getLBFGS(energy,commands,MINIMIZE);
  double               initialTemperature   = commands.firstWordFilter(MCM).secondWord(INITIAL_TEMPERATURE).thirdWordDouble();
  double               finalTemperature     = commands.firstWordFilter(MCM).secondWord(FINAL_TEMPERATURE).thirdWordDouble();
  int                  nSteps               = commands.firstWordFilter(MCM).secondWord(MAX_STEPS).thirdWordInt();
  TemperatureGenerator temperatureGenerator = new TemperatureGenerator(initialTemperature, finalTemperature,  nSteps);
  return new MCM(this, energy, minimizer, perturbation, temperatureGenerator, nSteps);
  }

}
