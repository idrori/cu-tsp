/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.molecularElements.atoms.Atom;
import meshi.parameters.ResidueType;

public class ChargedResidueSideChainFilter implements Filter {
    public static final ChargedResidueSideChainFilter filter = new ChargedResidueSideChainFilter();

    public boolean accept(Object obj) {
        Atom atom = (Atom) obj;
        return (!atom.type().backbone()) &&
                (atom.residue().type().equals(ResidueType.ASP) || atom.residue().type().equals(ResidueType.GLU) || atom.residue().type().equals(ResidueType.LYS));
    }
}