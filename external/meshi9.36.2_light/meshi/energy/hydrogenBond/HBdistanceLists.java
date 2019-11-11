package meshi.energy.hydrogenBond;

import meshi.geometry.*;
import meshi.util.filters.Filter;
/**
 **/
public class HBdistanceLists extends DistanceLists {
    final public Filter filter;
    public HBdistanceLists(Filter filter) {
        super(100);
        this.filter = filter;
    }

    private GoodResiduesForHB goodResiduesForHB = new GoodResiduesForHB();

}

