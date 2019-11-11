/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.util.filters.Filter;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.ResidueType;

/**
 *
 */
public class HydrophobicSideChains implements Filter {
    boolean includeCB;
    public HydrophobicSideChains(boolean includeCB) {
        this.includeCB = includeCB;
    }
    public boolean accept(Object o) {
        Atom atom = (Atom) o;
        if (atom.type().backboneC()) return false;
        if (atom.type().backboneCA()) return false;
        if ((atom.name().equals("CB")) && (!includeCB) )
            return false; // For compatibility with the backboneFilter

        ResidueType residueType = atom.residue().type;
        if (residueType == ResidueType.ARG) return false;
        if (residueType == ResidueType.CYS) return false;
        if (residueType == ResidueType.ASP) return false;
        if (residueType == ResidueType.GLU) return false;
        if (!includeCB) if (residueType == ResidueType.HIS) return false;
        if (residueType == ResidueType.LYS) return false;
        if (residueType == ResidueType.ASN) return false;
        if (residueType == ResidueType.PRO) return false;
        if (residueType == ResidueType.GLN) return false;
        if ((residueType == ResidueType.SER)) return false;
        if (atom.type().isCarbon()) return true;
        return false;
    }
}