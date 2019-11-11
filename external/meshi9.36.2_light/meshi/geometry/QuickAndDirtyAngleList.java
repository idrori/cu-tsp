/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

/**
 * A list of quick and dirty  angles. These angles may use approximate
 * angle calculations if ArcCos.useFastArcCos() is called.
 */
public class QuickAndDirtyAngleList extends AngleList {
    public QuickAndDirtyAngleList() {
        super();
    }

    public QuickAndDirtyAngleList(AtomPairList bonds, DistanceMatrix distanceMatrix) {
        super(bonds, distanceMatrix);
    }

    public Angle getAngle(AtomPair atomPair1, AtomPair atomPair2, DistanceMatrix distanceMatrix) {
        return new QuickAndDirtyAngle(atomPair1, atomPair2, distanceMatrix);
    }


    public static QuickAndDirtyAngleList getCaAnglesQuickAndDirty(Protein protein, DistanceMatrix distanceMatrix) {
        //AtomList CAs = protein.getAtoms().CAFilter();
        AtomList CAs = protein.atoms().CAnoImageFilter();
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
        return new QuickAndDirtyAngleList(CaPairs, distanceMatrix);
    }

    /**
     * Returns a sub-list containing angles that have a known name
     */
    public QuickAndDirtyAngleList namedFilterQD() {
        QuickAndDirtyAngleList out = new QuickAndDirtyAngleList();
        for (Angle angle : this)
            if (isNamed(angle)) out.add(angle);
        return out;
    }
}
