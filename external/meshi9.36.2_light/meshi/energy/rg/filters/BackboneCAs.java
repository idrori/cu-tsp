/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg.filters;

import meshi.util.filters.Filter;
import meshi.molecularElements.atoms.Atom;

/**
 * .
 */
public class BackboneCAs implements Filter {
    public boolean accept(Object o) {
        Atom atom = (Atom) o;
        if (atom.nowhere()) return false;
        if (atom.type().backboneCA()) return true;
        return false;
    }

}
