package meshi.util.filters;

import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;

public class CbAndGlyFilter implements Filter {
    public static final Filter filter = new CbAndGlyFilter();
    public boolean accept(Object obj) {
        if (CbFilter.filter.accept(obj)) return true;
        if (((Atom) obj).type() == AtomType.GCA) return true;
        return false;
    }
}
