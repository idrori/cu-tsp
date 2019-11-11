/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.util.filters.Filter;
import meshi.parameters.AtomType;
import meshi.molecularElements.atoms.Atom;

public class PolarSideChains implements Filter {

    public boolean accept(Object o) {
        Filter filter = new AtomType.ChargedFilter();
        Atom atom = (Atom) o;
        if (filter.accept(o)) return true;
        if (    (atom.type() == AtomType.NND) ||
                (atom.type() == AtomType.NOD) ||
                (atom.type() == AtomType.QNE) ||
                (atom.type() == AtomType.QOE) ||
                (atom.type() == AtomType.TOG) ||
                (atom.type() == AtomType.SOG)) return true;
        //
        return false;
    }
}
