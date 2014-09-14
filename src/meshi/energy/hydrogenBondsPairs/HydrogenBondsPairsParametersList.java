package meshi.energy.hydrogenBondsPairs;

import meshi.energy.Parameters;
import meshi.energy.*;
import meshi.util.filters.Filter;
import meshi.energy.simpleEnergyTerms.*;

public abstract class HydrogenBondsPairsParametersList  extends ParametersList{

    private final Parameters nonExsisting = new HydrogenBondsPairsParameters(1,1,1,1,1,1,-5); //todo check -5 (use to be 0)
    public int searchCounter =0;
    //private HydrogenBondsPairsParametersList coilParametersList;

    //---------------------------------- constructor -------------------------------
    public HydrogenBondsPairsParametersList (){
        super();                //create an empty list with filter IsParameter
    }



    public HydrogenBondsPairsParametersList(String parametersFileName) {
        super(parametersFileName ,true);  //true means it is sortable parameter list using compare
    }

    //---------------------------------- methods -----------------------------------

     public Parameters parameters(Object baseElement) {
        searchCounter ++;
        Parameters key = getKey(baseElement);
        Parameters params =   getParameters(key);
           if (params != null)
                return params;
         return nonExsisting ;
     }

    /**
     *
     * @param baseElement should be instance of PairOfHydrogenBondsElements
     * @return HydrogenBondsPairsParameters with 6 values to use un serching the freqency value
     */
    public Parameters getKey(Object baseElement) {
        //if (!(baseElement instanceof PairOfHydrogenBondsElements))
        //    throw new RuntimeException("problem with getKey in HBPEnergy"+ baseElement);
        PairOfHydrogenBondsElements element = (PairOfHydrogenBondsElements)baseElement;
        return new HydrogenBondsPairsParameters(element.seqDistance1,
                                                element.seqDistance2,
                                                element.hhSeqDistance,
                                                element.ooSeqDistance,
                                                element.hoSeqDistance,
                                                element.ohSeqDistance);
    }

    public Parameters createParameters(String line) {
        return new HydrogenBondsPairsParameters(line);
    }

    private static class NoElements implements Filter{
            public boolean accept(Object obj){
                    System.out.println("HydrogenBondsPairsParametersList:NoElements try to add an element");
                     return false;
            }
    }
}


