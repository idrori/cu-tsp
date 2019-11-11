package meshi.energy.rg.filters;

import meshi.molecularElements.atoms.Atom;
import meshi.util.filters.Filter;

/**
 * Created by chen on 22/11/2015.
 */
public class HereFilter implements Filter {
    public HereFilter() {};
    public boolean accept(Object obj) {
        return !((Atom) obj).nowhere();
    }
}
