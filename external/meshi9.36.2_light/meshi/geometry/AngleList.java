/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.util.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

import java.util.*;

/**
 **/
public class AngleList extends ArrayList<Angle> implements Updateable {
    private int numberOfUpdates = 0;

    public AngleList() {
        super();
    }

    public AngleList(AtomPairList bonds, DistanceMatrix distanceMatrix) {
        super();
        AtomPair bond1, bond2;
        for (int iBond1 = 0; iBond1 < bonds.size(); iBond1++) {
            bond1 = bonds.get(iBond1);
            if (distanceMatrix.distance(bond1) != null) {
                for (int iBond2 = 0; iBond2 < iBond1; iBond2++) {
                    bond2 = bonds.get(iBond2);
                    if (bond2.atom1().nowhere() | bond2.atom2().nowhere())
                        continue;
                    if (distanceMatrix.distance(bond2) != null) {
                        if (bond1.sharedAtom(bond2) != null)
                            add(getAngle(bond1, bond2, distanceMatrix));
                    }
                }
            }
        }
    }
   public void update(int numberOfUpdates) throws UpdateableException {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            int size = size();
            for (int i = 0; i < size; i++) {
                angleAt(i).update(numberOfUpdates);
            }
            this.numberOfUpdates++;
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with AngleList.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }

    public Angle angleAt(int i) {
        return (Angle) get(i);
    }

    /**
     * Returns a sub-list containing angles that have a known name     *
     */
    public AngleList namedFilter() {
        Iterator angles = iterator();
        AngleList out = new AngleList();
        for (Angle angle : this)
            if (isNamed(angle)) out.add(angle);
        return out;
    }

    public boolean isNamed(Angle angle) {
        if (angle.getAngleName().compareTo("") != 0)
            return true;
        else
            return false;
    }


    public AtomList atomList() {
        if (size() < 1) throw new RuntimeException("Cannot create atom list for an empty angle list.");
        AtomList list = new AtomList(get(0).atom1.molecularSystem);
        Atom atom1, atom2, atom3;
        for (Angle angle : this) {
            atom1 = angle.atom1();
            atom2 = angle.atom2();
            atom3 = angle.atom3();
            if (!list.contains(atom1)) {
                list.add(atom1);
            }
            if (!list.contains(atom2)) {
                list.add(atom2);
            }
            if (!list.contains(atom3)) {
                list.add(atom3);
            }
        }
        return list;
    }


    public static AngleList getCaAngles(Protein protein, DistanceMatrix distanceMatrix) {
        AtomList CAs = protein.atoms().CAFilter();
        AtomPairList CaPairs = new AtomPairList();

        Atom atom1, atom2;
        for (int iAtom1 = 0; iAtom1 < CAs.size(); iAtom1++) {
            atom1 = CAs.get(iAtom1);
            for (int iAtom2 = 0; iAtom2 < iAtom1; iAtom2++) {
                atom2 = CAs.get(iAtom2);
                if (Math.abs(atom1.residueNumber() - atom2.residueNumber()) < 2) {
                    if (atom1.residueNumber() < atom2.residueNumber())
                        CaPairs.add(new AtomPair(atom1, atom2));
                    else
                        CaPairs.add(new AtomPair(atom2, atom1));
                }
            }

        }
        return new AngleList(CaPairs, distanceMatrix);
    }


    public Angle getAngle(AtomPair atomPair1, AtomPair atomPair2, DistanceMatrix distanceMatrix) {
        return new Angle(atomPair1, atomPair2, distanceMatrix);
    }


    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
        for (Angle angle : this)
            angle.setNumberOfUpdates(numberOfUpdates);
    }
}

