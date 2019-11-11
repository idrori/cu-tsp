/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB;

import meshi.molecularElements.atoms.Atom;
import meshi.util.filters.Filter;

public class OneOfTwenty implements Filter {
    public boolean accept(Object obj) {
        if (!(obj instanceof Atom))
            throw new RuntimeException("Weird parameter to OneOfTwenty.accept " + obj);
        Atom atom = (Atom) obj;
        if (atom.residueName().equals("ALA")) return true;
        if (atom.residueName().equals("CYS")) return true;
        if (atom.residueName().equals("ASP")) return true;
        if (atom.residueName().equals("GLU")) return true;
        if (atom.residueName().equals("PHE")) return true;
        if (atom.residueName().equals("GLY")) return true;
        if (atom.residueName().equals("HIS")) return true;
        if (atom.residueName().equals("ILE")) return true;
        if (atom.residueName().equals("LYS")) return true;
        if (atom.residueName().equals("LEU")) return true;
        if (atom.residueName().equals("MET")) return true;
        if (atom.residueName().equals("ASN")) return true;
        if (atom.residueName().equals("PRO")) return true;
        if (atom.residueName().equals("GLN")) return true;
        if (atom.residueName().equals("ARG")) return true;
        if (atom.residueName().equals("SER")) return true;
        if (atom.residueName().equals("THR")) return true;
        if (atom.residueName().equals("VAL")) return true;
        if (atom.residueName().equals("TRP")) return true;
        if (atom.residueName().equals("TYR")) return true;
        return false;
    }
}
	    
