/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg.filters;

import meshi.molecularElements.atoms.Atom;
import meshi.util.filters.Filter;

/**
 *
 */
public class HeavyAtoms implements Filter {
    public boolean accept(Object o) {
        Atom atom = (Atom) o;
        return ((!atom.isHydrogen()) && (!atom.nowhere()));
    }

}

