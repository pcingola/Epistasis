package meshi.symmetryComplex.topologyMap;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.Residue;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;
import meshi.symmetryComplex.utils.GJFilters;

/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Dec 23, 2008
 * Time: 2:45:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class Topology implements KeyWords {

    /*public static final int[]
N_TER = { 1, 18 },
M1 = { 19, 40 },
E1 = { 41, 75 },
M2 = { 76, 96 },
IN = { 97, 130 },
M3 = { 131, 152 },
E2 = { 153, 188 },
M4 = { 189, 209 },
C_TER = { 210, 999 };*/
    //private Part predictedPart, loopPart, absentPart;
    private Loop [] loops;
    private final String LOOP = "loop";
          //
    public Topology(CommandList commands) {

              //predictedPart = new Part(BoundariesMap.predictedMap(), residues);
       //       absentPart = new Part(BoundariesMap.absentMap(), residues);
         //     loopPart = new Part(BoundariesMap.loopMap(), residues);
//more about loops
              loops = new Loop[BoundariesMap.loopMap().length];
              for (int i = 0; i < BoundariesMap.loopMap().length; i++)
                     loops[i] = new Loop(commands.firstWordFilter(LOOP+String.valueOf(i+1)),i+1);
    }

    //public Part predictedPart(){return predictedPart;}
    //public Part absentPart(){return absentPart;}
    //public Part loopPart(){return loopPart;}

    public int [][] predictedMap() {return BoundariesMap.frozenMap();}
    public int [][] absentMap() {return BoundariesMap.absentMap();}
    public int [][] loopMap() {return BoundariesMap.loopMap();}

    //public ResidueList predictedPartResidues(){return predictedPart.residues();}
    //public ResidueList absentPartResidues(){return absentPart.residues();}
    //public ResidueList loopsResidues(){return loopPart.residues();}

    public Loop [] loops(){return loops;}
    public Loop loop(int i){return loops[i-1];}
    public int [] loopBounds(int i){return BoundariesMap.loopMap(i);}

/*
  public ResidueList residuesOfPredictedPart(ResidueList residues) {
      return residues.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.frozenMap()));
  }

  public ResidueList residuesOfAbsentPart(ResidueList residues) {
      return residues.filter(new GJFilters.ResidueRangeFilter(BoundariesMap.absentMap()));
  }

  public ResidueList residuesOfPredictedPartChainA(SymmetricComplex gj) {
        return gj.getSource().filter(new GJFilters.ResidueRangeFilter(BoundariesMap.frozenMap()));
  }

  public ResidueList residuesOfPredictedPartChainNum(SymmetricComplex gj, int chainNum) {
          return gj.chain(chainNum).filter(new GJFilters.ResidueRangeFilter(BoundariesMap.frozenMap()));
  }


    public ResidueList residuesOfAbsentPartChainA(SymmetricComplex gj) {
          return gj.getSource().filter(new GJFilters.ResidueRangeFilter(BoundariesMap.absentMap()));
    }

    public ResidueList residuesOfAbsentPartChainNum(SymmetricComplex gj, int chainNum) {
            return gj.chain(chainNum).filter(new GJFilters.ResidueRangeFilter(BoundariesMap.absentMap()));
    }

    public ResidueList residuesOfLoopsChainA(SymmetricComplex gj) {
          return gj.getSource().filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap()));
    }

    public ResidueList residuesOfLoopsChainNum(SymmetricComplex gj, int chainNum) {
            return gj.chain(chainNum).filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap()));
    }

    public ResidueList residuesOfLoopNumChainNum(SymmetricComplex gj, int loopNum, int chainNum) {
            return gj.chain(chainNum).filter(new GJFilters.ResidueRangeFilter(BoundariesMap.loopMap(loopNum)));
    }
 */

    /**
     * Returns '1' if the residue number belongs to loop E1, '2' if it belongs
     * to loop E2 and '0' otherwise.
     */
    public int whichLoop(final int residueNumber) {
        for (int i=1; i<=loops.length; i++)
        if (residueNumber >= loopBounds(i)[0] & residueNumber <= loopBounds(i)[1])
            return i;
        return 0;
    }

    public boolean inLoop(final Residue residue) {
        if (residue == null)
            return false;
        int residueNumber = residue.number();
        for (int i=1; i<=loops.length; i++)
          if (residueNumber >= loopBounds(i)[0] & residueNumber <= loopBounds(i)[1])
            return true;
        return false;
    }

     public boolean inLoop(final int residueNumber) {
         for (int i=1; i<=loops.length; i++)
           if (residueNumber >= loopBounds(i)[0] & residueNumber <= loopBounds(i)[1])
             return true;
         return false;
     }

     public int loopFirstAtomNumber(int loopNumber) {
        if (loopNumber < 1 || loopNumber >= loops.length+1)
//            return -1;
            throw new RuntimeException("\nTopology.loopFirstAtomNumber tries to call to unexisted loop\n");
        return loopBounds(loopNumber)[0];

    }

     public int loopLastAtomNumber(int loopNumber) {
         if (loopNumber < 1 || loopNumber >= loops.length+1)
//            return -1;
             throw new RuntimeException("\nTopology.loopFirstAtomNumber tries to call to unexisted loop\n");
         return loopBounds(loopNumber)[1];
    }
}
