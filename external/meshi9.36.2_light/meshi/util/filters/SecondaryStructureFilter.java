/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.molecularElements.atoms.Atom;
import meshi.parameters.SecondaryStructure;

/**
 *
 */
public class SecondaryStructureFilter implements Filter {
    public boolean accept(Object o) {
        Atom atom = (Atom) o;
        SecondaryStructure secondaryStructure = atom.residue().getSecondaryStructure();
        if ((secondaryStructure == SecondaryStructure.HELIX) || (secondaryStructure == SecondaryStructure.SHEET))
            return true;
        return false;
    }
}
