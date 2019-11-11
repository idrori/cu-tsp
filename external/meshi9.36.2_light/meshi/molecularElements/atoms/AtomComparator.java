/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.atoms;

import meshi.molecularElements.*;

import java.util.*;

public class AtomComparator implements Comparator {
    public int compare(Object obj1, Object obj2) {
        Atom atom1 = (Atom) obj1;
        Atom atom2 = (Atom) obj2;
        if (atom1.residue().number() > atom2.residue().number()) return 1;
        if (atom1.residue().number() < atom2.residue().number()) return -1;
        return (atom1.name().compareTo(atom2.name()));
    }
}