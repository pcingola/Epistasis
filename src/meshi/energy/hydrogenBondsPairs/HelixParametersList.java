package meshi.energy.hydrogenBondsPairs;

import meshi.util.filters.Filter;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 06/02/2006
 * Time: 12:34:10
 *  A specific HydrogenBondsPairsParametersList for helices
 */
public class HelixParametersList extends HydrogenBondsPairsParametersList {
       //------------------------------------ constructors ----------------------------
      public HelixParametersList(){
          super();
     }

     public HelixParametersList(String parametersFileName) {
        super(parametersFileName);
    }

    //---------------------------------------------------------------------------------------------
    private static class IsHydrogenBondsPairsParameters implements Filter {
    public boolean accept(Object obj) {
        return (obj instanceof HydrogenBondsPairsParameters);
    }
    }
}
