/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology;

import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.filters.Filter;

import java.util.Iterator;

public class NonFrozensRadiusFilter implements Filter {

    AtomList /*getAtoms,*/nonFrozens;
    double radius;
    int bondDepth;

    /*public NonFrozensRadiusFilter(AtomList getAtoms, double radius, int bondDepth) {

        *//*this.getAtoms = getAtoms;*//*
        this.radius = radius;
        this.bondDepth = bondDepth;

        nonFrozens = new AtomList();
        Atom atom;
        Iterator atomListIter = getAtoms.iterator();
        while ((atom = (Atom) atomListIter.next()) != null) {
            if (!atom.frozen()) nonFrozens.add(atom);
        }
    }*/

    public boolean accept(Object obj) {
        Atom atom = (Atom) obj;
        /*if (!getAtoms.contains(atom))
        throw new RuntimeException("atom "+atom+" is not in the filter list");
*/
        if (!atom.frozen()) return true;

        if (isConnectedToNonFrozen(atom, bondDepth)) return true;

        if (inRadiusFromNonFrozen(atom)) return true;

        return false;
    }

    private boolean isConnectedToNonFrozen(Atom atom, int depth) {

        if (!atom.frozen()) return true;
        if (depth == 0) return false;

        Atom bonded;
        Iterator bondeds = atom.bonded().iterator();
        while ((bonded = (Atom) bondeds.next()) != null) {
            if (isConnectedToNonFrozen(bonded, depth - 1)) return true;
        }

        return false;
    }

    private boolean inRadiusFromNonFrozen(Atom atom) {

        Atom nonFrozen;
        Iterator nonFrozensIter = nonFrozens.iterator();
        while ((nonFrozen = (Atom) nonFrozensIter.next()) != null) {
            if (atom.distanceFrom(nonFrozen) <= radius) return true;
        }
        return false;
    }
}
