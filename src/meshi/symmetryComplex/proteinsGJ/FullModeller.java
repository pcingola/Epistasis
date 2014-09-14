package meshi.symmetryComplex.proteinsGJ;

import java.io.IOException;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.CooperativeAtomicPairwisePMFSummaCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativePropensityCreator;
import meshi.energy.simpleEnergyTerms.inflate.simpleInflate.SimpleInflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.simpleInflate.SimpleInflate;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ResidueCreator;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.BackboneResidue;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.ca.CaResidue;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.Command;
import meshi.util.Scmod;
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

public class FullModeller extends BackboneCompleterAndMinimizer {


    public FullModeller(String pdbTemplate, String sequence, CommandList commands,
                        String outputFileMask) throws IOException {

        super(sequence, commands, outputFileMask);
// to determine if pdbTemplate the file of CA for TM helices or this is backbone template for full GJ
        AtomList tmp = new AtomList(pdbTemplate);
        boolean isCaStage = true;
        if ( !(tmp.get(0).isCarbon() && tmp.get(0).isCarbon())) isCaStage = false;


    if (isCaStage) { // if only CA positions of TM helixes are known
        setGJ(pdbTemplate);

        if (!toContinue("CA")) return;
        caRelax();
        if (!toContinue("CA")) return;
        
        restoreAll(creator);
        generateFiles(gj,creator);
         timeStamp();
        generateLoopData("CA");
       if (!toContinue("CA")) return;

        backboneMinimization();
    }
    else {
        backboneMinimization(pdbTemplate);
    }

     //   update();
        generateFiles(gj,creator);
        timeStamp();

        // complete sidechains?
        creator = ResidueExtendedAtomsCreator.creator;
        new MolecularSystem();
//        Protein schp = new Protein(gj.getSource().atoms(), creator);
//        Protein schp = new Protein(gj.getSource().atoms(), new ResidueExtendedAtoms(ADD_ATOMS));
//to freeze
        //schp.atoms().freeze(GJFilters.frozenResiduesFilter);

	Protein schp = new Protein(gj.name());
	ResidueExtendedAtomsCreator creatorAll = new ResidueExtendedAtomsCreator();
	new MolecularSystem();
	Chain backboneChain = gj.chains().get(0);
	Chain chain = new Chain(backboneChain.name(),schp);
        for (Object aBackboneChain : backboneChain) {
            Residue residue = (Residue) aBackboneChain;
            if (residue.dummy()) chain.add(residue);
            else {
                Residue newResidue = creatorAll.create(residue.atoms(), residue.ID(), residue.mode);
                newResidue.setSecondaryStructure(residue.secondaryStructure());
                chain.add(newResidue);
            }
        }
        schp.addChain(chain);
        Scmod.scmod(commands, schp, 1,50);
        
        new MolecularSystem();
        SymmetricComplex gjAll = new SymmetricComplex(schp.chains().get(0).atoms(),creator,super.alignment,getTransformations(RELATIVE_ORIENTATION,numberOfChains));
        generateFiles(gjAll,creator);
        timeStamp();

        gj = makeActualSymmetryComplex(gjAll, creator);

         // minimize full model?
        fullMinimization();

        //update();
        generateFiles(gj,creator);
        timeStamp();
        generateData("All");
        generateLoopData("All");
    }


    protected void fullMinimization() throws IOException {
        System.out.println("MCM");
        gj.printChainsNames();

//            SolvateCBCreator solvateCBCreator = new SolvateCBCreator();
        TetherCreator tetherCreator = new TetherCreator();
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
                summaCreator,
                new CooperativeAtomicPairwisePMFSummaCreator(summaCreator),
                new SolvateCreatorHBforMinimization(1.0 , 1.0, 1.0, 1.0,commands),
                hydrogenBondsCreator,
                new HydrogenBondsPairsCreator(hydrogenBondsCreator),
                new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
                new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
                //tetherCreator,

               new EDMEnergyCreator(GJFilters.chainA_Ca_loops_Atom_Filter),
     //         getCysCBetasDistConstraintsCreator(),
               new CylinderCreator(commands, "innerBarrel", new GJFilters.AndFilter(new CaAtomFilter(), innerBarrelFilter)),
               new CylinderCreator(commands, "outerBarrel", new GJFilters.AndFilter(new CaAtomFilter(), outerBarrelFilter)),
               inflateCreator,
            };
//            DistanceMatrixNoImageAtoms distanceMatrix =toDo Solvate
                    //new DistanceMatrixNoImageAtoms(gj, solvateCBCreator.getParametersList(commands));
        DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj.atoms().molecularSystem(),
                            5.5, 2.0, DistanceMatrix.DEFAULT_EDGE(5.5,2.0), 4, gj.imageResidues());
            energy = new TotalEnergy(gj, distanceMatrix, creators, commands);

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
    Perturbation         perturbation         = new Perturbation(energy,commands, inflate, model,0);
    Minimizer            minimizer            = Utils.getLBFGS(energy,commands,MINIMIZE);
    double               initialTemperature   = commands.firstWordFilter(MCM).secondWord(INITIAL_TEMPERATURE).thirdWordDouble();
    double               finalTemperature     = commands.firstWordFilter(MCM).secondWord(FINAL_TEMPERATURE).thirdWordDouble();
    int                  nSteps               = commands.firstWordFilter(MCM).secondWord(MAX_STEPS).thirdWordInt();
    TemperatureGenerator temperatureGenerator = new TemperatureGenerator(initialTemperature, finalTemperature,  nSteps);
    return new MCM(this, energy, minimizer, perturbation, temperatureGenerator, nSteps);
    }

}
