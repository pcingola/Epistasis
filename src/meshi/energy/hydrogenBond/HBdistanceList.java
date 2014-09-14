package meshi.energy.hydrogenBond;

import meshi.geometry.*;
import meshi.util.filters.Filter;
/**
 **/
public class HBdistanceList extends DistanceList {

    public HBdistanceList(Filter fi) {
        super(100,fi);
    }

    private GoodResiduesForHB goodResiduesForHB = new GoodResiduesForHB();

}

