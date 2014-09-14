package meshi.symmetryComplex.proteinsGJ;

import static meshi.symmetryComplex.utils.GJUtils.*;
import static meshi.symmetryComplex.utils.LoopUtils.*;

import meshi.energy.TotalEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.alphaTorsion.AlphaTorsionCreator;
import meshi.energy.simpleEnergyTerms.alphaAngle.AlphaAngleCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativePropensityCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.CooperativeAtomicPairwisePMFSummaCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HBondList;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.energy.hydrogenBond.GoodResiduesForHB;
import meshi.util.*;
import meshi.util.filters.Filter;
import meshi.util.file.MeshiWriter;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.ca.CaResidue;
import meshi.sequences.Sequence;
import meshi.sequences.SequenceAlignment;
import meshi.applications.prediction.homology.CaAtomFilter;
import meshi.applications.prediction.PredictionSequenceAlignment;
import meshi.optimizers.LBFGS;
import meshi.optimizers.Optimizer;
import meshi.geometry.Distance;
import meshi.symmetryComplex.geometryImage.DistanceMatrixNoImageAtoms;
import meshi.geometry.DistanceMatrix;
import meshi.symmetryComplex.topologyMap.Topology;
import meshi.symmetryComplex.topologyMap.BoundariesMap;
import meshi.symmetryComplex.utils.GJFilters;
import meshi.symmetryComplex.utils.LoopUtils;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplexCompleter;
import meshi.symmetryComplex.energy.cylinder.CylinderCreator;
import meshi.symmetryComplex.energy.edmEnergy.EDMEnergyCreator;
import meshi.symmetryComplex.energy.SymmetryCreator;
import meshi.symmetryComplex.transformations.Transformation;

import java.io.IOException;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.GregorianCalendar;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Dec 23, 2008
 * Time: 2:36:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class GJModeller implements KeyWords{
    public SymmetricComplex gj;
    protected TotalEnergy energy;
    protected CommandList commands;
    protected char stage = 'a'; // The files generated in stage 'a' are practically useless, but are important because of random numbers generation.
    protected String outputFileMask;
    protected Sequence sequence;
    public Protein caHelicesProtein;
    protected static Topology topology;
    protected static final double DISTANCE_AS_NEIGHBOUR_CHAINS = 10;
    protected static int numberOfChains;
    protected SequenceAlignment alignment;
    protected ResidueCreator creator;
    public ResidueCreator creator(){return creator;}


    public GJModeller(String sequence, CommandList commands, String outputFileMask) throws IOException {
        numberOfChains = BoundariesMap.numberOfChains;
        if (numberOfChains < 1)
        throw new RuntimeException("Something weird with BoundariesMap.numberOfChains or with BoundariesMap initialization");
        this.commands = commands;
        this.outputFileMask = outputFileMask;
        this.sequence = new Sequence(sequence, "1txhA");
        alignment = new PredictionSequenceAlignment(commands, SECONDARY_STRUCTURE, this.sequence);
        topology = new Topology(commands);
    }


    public GJModeller(String tmCasPdbLocation, String sequence, CommandList commands,
                      String outputFileMask) throws IOException {
        this(sequence, commands, outputFileMask);
        setGJ(tmCasPdbLocation);
    }

    protected void setGJ(String tmCasPdbLocation) throws IOException{
        creator = new CaResidue("creator");
        this.caHelicesProtein = new Protein(new AtomList(tmCasPdbLocation).splitToChains()[0], creator);

        new MolecularSystem();
         //caProtein = the Protein without absent atoms, but with NOWHERE atoms
        Protein caProtein = new Protein(sequence,"1thx","A",creator);
        positionTMCas(caHelicesProtein, caProtein);

     //   new MolecularSystem();
        SymmetricComplexCompleter gjCompleter = new SymmetricComplexCompleter(caProtein.chain().atoms(), creator,
                getTransformations(RELATIVE_ORIENTATION,numberOfChains),commands, topology);

        gj = makeActualSymmetryComplex(gjCompleter,creator);

        generateFiles(gj, creator);
        timeStamp();

    }

    protected SymmetricComplex makeActualSymmetryComplex(SymmetricComplex gjAll, ResidueCreator creator) throws IOException {
     /*
        new MolecularSystem();
        gj = new SymmetricComplex(gjCompleter.getSource().atoms(), creator, alignment,
                getTransformations(RELATIVE_ORIENTATION,numberOfChains));

        generateFiles(gj, creator);
        timeStamp();
        */
        boolean[] chainsBooleans = whichChainsInZone(DISTANCE_AS_NEIGHBOUR_CHAINS, gjAll.atoms());
        //getChainLetters(chainsBooleans)
        Transformation[] partialTransformations = getPartialTransformations(RELATIVE_ORIENTATION, numberOfChains, chainsBooleans);

       return  new SymmetricComplex(gjAll.getSource().atoms(), creator, alignment,
                partialTransformations, getChainLetters(chainsBooleans));
        //*/
    }



    protected void restoreAll(ResidueCreator creator){
        if (gj.chains().size() == BoundariesMap.numberOfChains) return;
        new MolecularSystem();
         gj = new SymmetricComplex(gj.getSource().atoms(), creator, alignment,
                      getTransformations(RELATIVE_ORIENTATION, numberOfChains));


    }

    protected boolean toContinue(String stage){                                 //toDo
         //if need to continue calculation
         System.out.println("\n#"+stage+"# Crosses of Loops");
         int numOfCrosses = 0;
  //           for (int ch1 = 0; ch1 < charLetters.length; ch1++)
         int num;
         //crosses for chain A
             for (int ch = 1; ch < gj.chainLetters().length(); ch++){

                 if (numberOfTwists(gj.atoms(),1) !=0 ) {
                     System.out.println("Twists in Loop E1");
                     return false;
                 }
                 if (numberOfTwists(gj.atoms(),2) !=0 ) {
                     System.out.println("Twists in Loop E2");
                     return false;
                 }

                num = LoopUtils.numberOfChainCrosses(gj.atoms(),"A",""+gj.chainLetters().charAt(ch));
                if (num != 0)
                    System.out.print(num+"(A,"+gj.chainLetters().charAt(ch)+")\t");
                numOfCrosses +=num;
             }
         System.out.println("Total number of crosses is "+numOfCrosses);

         return (numOfCrosses == 0);
     }

    protected void caRelax() throws IOException {
        System.out.println("C-Alpha relaxation");
        gj.printChainsNames();

      //  gj.updateBondsList();
        EnergyCreator[] creators = {
                new SymmetryCreator(),
                new BondCreator(),
                new AlphaAngleCreator(),
                new PlaneCreator(),
                new AlphaTorsionCreator(),
                new ExcludedVolCreator(),                
    //            getCysCAlphasDistConstraintsCreator(),
                new CylinderCreator(commands, "innerBarrel", innerBarrelFilter),
                new CylinderCreator(commands, "outerBarrel", outerBarrelFilter),

        };

        DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj.atoms().molecularSystem(),
                            4.5, 1.0, DistanceMatrix.DEFAULT_EDGE(4.5,1.0),4,null);
        //DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj.molecularSystem(), gj.fullMolecularSystem(),
          //                  4.5, 1.0, DistanceMatrix.DEFAULT_EDGE(4.5,1.0),4,null);

        energy = new TotalEnergy(gj, distanceMatrix, creators, commands);
        LBFGS lbfgs = Utils.getLBFGS(energy, commands,RELAX);
        try {
            Optimizer.OptimizerStatus os = lbfgs.run();
            System.out.println("Optimization "+os);
        }
        catch (Exception ex) {
            energy.test();
            throw new RuntimeException(ex);
        }


    }


    protected void backboneRelax() throws IOException {
        System.out.println("Backbone relaxation");
        gj.printChainsNames();

   //     SolvateCBCreator solvateCBCreator =  new SolvateCBCreator();
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

//                solvateCBCreator,
/*
                new ExcludedVolCreator(),
                new TwoTorsionsCreator(),
 //                new SolvateCreator(1.0, 1.0),
        //        new FastSolvateCreator(0.001, 1.0),
                hbc,
                new HydrogenBondsPairsCreator(hbc),
                new HBondsPunishOHNAngleCreator(hbc),
                new HbondsPunishHOCAngleCreator(hbc),
        //        new HydrogenBondsPlaneCreator(),
  */
                ramachCreator,
                new CooperativeRamachandranCreator(ramachCreator),
//                new CompositePropensityCreator(),
                propensityCreator,
                new CooperativePropensityCreator(propensityCreator),
                //new ExcludedVolCreator(),
                summaCreator,
               new CooperativeAtomicPairwisePMFSummaCreator(summaCreator),
                hydrogenBondsCreator,
                new HydrogenBondsPairsCreator(hydrogenBondsCreator),
                new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
                new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
                //new SolvateCreatorHBforMinimization(1.0 , 1.0, 1.0, 1.0,commands),  //All
//                new SimpleInflateCreator(),
                //new TetherCreator(),
                new EDMEnergyCreator(GJFilters.chainA_Ca_loops_Atom_Filter),
//                getCysCBetasDistConstraintsCreator(),
                new CylinderCreator(commands, "innerBarrel", new GJFilters.AndFilter(new CaAtomFilter(), innerBarrelFilter)),
                new CylinderCreator(commands, "outerBarrel", new GJFilters.AndFilter(new CaAtomFilter(), outerBarrelFilter)),

        };
//        DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj, solvateCBCreator.getParametersList(commands));
        DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj.atoms().molecularSystem(),
                                  5.5, 1.0, DistanceMatrix.DEFAULT_EDGE(5.5,1.0), 4, gj.imageResidues());
        //DistanceMatrixNoImageAtoms distanceMatrix = new DistanceMatrixNoImageAtoms(gj.molecularSystem(), gj.fullMolecularSystem(),
          //                        5.5, 1.0, DistanceMatrix.DEFAULT_EDGE(5.5,1.0), 4, gj.actualImageResidues());

        energy = new TotalEnergy(gj, distanceMatrix, creators, commands);

//        distanceMatrix.testNonBondedList();
      //  energy.test();
  //      energy.findCriminalEnergyTerm(gj.atoms().atomAt(754));
    //    energy.findCriminalEnergyTerm(gj.atoms().atomAt(754));
      //  energy.findCriminalEnergyTerm(gj.atoms().atomAt(754));

   LBFGS lbfgs =  Utils.getLBFGS(energy, commands,RELAX);
        try {
            Optimizer.OptimizerStatus os = lbfgs.run();
            System.out.println("Optimization "+os);
        }
        catch (Exception ex) {
            energy.test();
            throw new RuntimeException(ex);
        }
    }


    public void generateData(String stage) throws IOException {
        MeshiWriter writer = new MeshiWriter(outputFileMask+"-"+stage+".dat");

        // Getting and writing the values of energy total and terms, as they appear
        // in a report.
        String report = energy.report(0);
        int i = report.indexOf('\n');
        if (i>=0) report = report.substring(i+1);
        StringTokenizer reportPieces = new StringTokenizer(report);
        reportPieces.nextToken(); // Dump the meaningless "iteration number".
        String token = reportPieces.nextToken();
        writer.print("EnergyReport\t"+token);
        System.out.print("EnergyReport"+stage+"\t"+token);

        while (reportPieces.hasMoreTokens())  {
            token = reportPieces.nextToken();
            writer.print("\t"+ token);
            System.out.print("\t"+token);
        }
        System.out.println();

        writer.print(getAverageDistanceBetweenCysCBconnexin("CB",gj.atoms().splitToChains()[0]));
//        writer.print("\nDistance Between Loops "+getAverageDistanceBetweenLoops(gj.atoms().splitToChains()[0]));

        // Getting and writing the number of H-bonds between chains and in total.
        boolean found = false;
        Iterator terms = energy.energyTerms().iterator();
        AbstractEnergy ae = null;
        while (!found) {
            ae = (AbstractEnergy) terms.next();
            if (ae instanceof HydrogenBondsEnergy)
                found = true;
        }
        HBondList hBonds = ((HydrogenBondsEnergy) ae).hBondList();
        writer.println("\t"+hBonds.size());
        writer.println();
        writer.close();
}


    public void generateLoopData(String stage) throws IOException {
          getDistancesBetweenCysCBofConnexin("CA",gj.atoms().splitToChains()[0]);
          getDistancesBetweenCysCBofConnexons("CA",
                          gj.atoms().splitToChains()[0],gj.atoms().splitToChains()[1],
                          gj.atoms().splitToChains()[6],gj.atoms().splitToChains()[7]);

        MeshiWriter writer = new MeshiWriter(outputFileMask+"-"+stage+".loop");
        int numOfParts = 5;
        int numTwists;
        writer.println("#Seed  Loop NumberOfTwists AverageAngle AnglesForEachPart ... AverageAngleWithoutFirst/LastPartOfLoop");
        for (int i=1; i<=BoundariesMap.loopMap().length; i++){
            numTwists = numberOfTwists(gj.atoms(),i);
            writer.print(i+"\t"+numTwists+"\t"+
                getAverageAngleForLoop(gj.atoms(),i));
            double [] partAngle = getAverageAngleForLoopParts(gj.atoms(),i, numOfParts);
            for (double angle : partAngle)
                    writer.print("\t"+angle);

            int neglect = 1;
            if (numOfParts > 5) neglect = 2;
            double trustSum = 0;
            for (int k = neglect; k < partAngle.length-neglect; k++)
                    trustSum += partAngle[k];
                    writer.print("\t"+trustSum/(partAngle.length-2*neglect));
            writer.println();
        }

        writer.println("\n#Crosses of Loops");
        int numOfCrosses = 0;
 //           for (int ch1 = 0; ch1 < charLetters.length; ch1++)
        int num;
        //crosses for chain A
            for (int ch = 1; ch < gj.chainLetters().length(); ch++){
               num = LoopUtils.numberOfChainCrosses(gj.atoms(),"A",""+gj.chainLetters().charAt(ch));
               if (num != 0)
                   writer.print(num+"(A,"+gj.chainLetters().charAt(ch)+")\t");
               numOfCrosses +=num;
            }
        writer.println(numOfCrosses);
        writer.close();
    }

    /**
     * Used by Window.
     */
    public final long startTime = new GregorianCalendar().getTimeInMillis();
    public void timeStamp() {
        long time = (new GregorianCalendar().getTimeInMillis() - startTime) / 1000;
        long hours = time / 3600, minutes = (time % 3600)/60, seconds = time % 60;
        System.out.println("timeStamp "+fixZero(hours)+":"+fixZero(minutes)+":"+fixZero(seconds));
    }

    private static String fixZero(long num) {
        return (num<10 ? "0" : "") + num;
    }

    public void generateFiles(SymmetricComplex sc, ResidueCreator creator) throws IOException {
        String filenamePrefix = outputFileMask+"-"+stage+"-";
        MeshiWriter pdb = new MeshiWriter(filenamePrefix+"atoms.pdb");
        AtomList atomList;
                      //  energy.;

        if (sc.chains().size() != BoundariesMap.numberOfChains){
            MolecularSystem current = MolecularSystem.currentMolecularSystem();

            new MolecularSystem();
      //      MolecularSystem newms = new MolecularSystem();
            SymmetricComplex s = new SymmetricComplex(gj.getSource().atoms(), creator,
                              getTransformations(RELATIVE_ORIENTATION, numberOfChains));
            atomList = s.atoms();
            MolecularSystem.setCurrentMolecularSystem(current);
            DistanceMatrix.terminator.reset(); 
        }
      else
          atomList = sc.atoms();

        atomList.print(pdb);
        pdb.close();
        System.out.println("\n*** Generated files, stage "+stage+" ***");
        stage++;
    }

    public String getChainLetters(boolean [] chainsBooleans){
        String chainLetters = "A";
        for (int i = 0; i < chainsBooleans.length; i ++)
            if (chainsBooleans[i]){
                String tmp = ""+(char)('B'+i);
                chainLetters += tmp;
            }
        return chainLetters;
    }


//--------------------------------filters for HB----------------------------------------------------
    public static class HBRangesFilter implements Filter {
        private Filter loop1Filter = new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(1));
        private Filter loop2Filter = new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(2));

        public boolean accept(Object obj) {
            Distance distance = (Distance) obj;
            Atom atom1 = distance.atom1(), atom2 = distance.atom2();
            return ((loop1Filter.accept(atom1) & loop1Filter.accept(atom1)) ^
                (loop2Filter.accept(atom1) & loop2Filter.accept(atom2)));
        }
    }

    public static class HBRangesAndGoodResiduesFilter extends GJFilters.AndFilter {
        public HBRangesAndGoodResiduesFilter() {
            super(new HBRangesFilter(), new GoodResiduesForHB(), new SpecialHBPairs());
        }
    }

    private static class SpecialHBPairs implements Filter {
        public boolean accept(Object obj) {
            Distance dis = (Distance) obj;
            return dis.atom1().core.status().image() ^ dis.atom2().core.status().image();
        }
    }


//--------------------------------filters for loops---------------------------------------

    public static Filter innerBarrelFilter = new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(2));
    public static Filter outerBarrelFilter = new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(1));



}

