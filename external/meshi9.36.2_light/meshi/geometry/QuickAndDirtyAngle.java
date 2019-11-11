/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

/**
 * Allows the use of approximated arc-cos function.
 * The use of the standard Math.acos function is the default.
 */

public class QuickAndDirtyAngle extends Angle {
    //private final ArcCos ARC_COS = new ArcCos();

    public QuickAndDirtyAngle(AtomPair atomPair1, AtomPair atomPair2,
                              DistanceMatrix distanceMatrix) {
        super(atomPair1, atomPair2, distanceMatrix);
    }

    public QuickAndDirtyAngle(Atom atom1, Atom atom2, Atom atom3, DistanceMatrix distanceMatrix) {
        this(new AtomPair(atom1, atom2), new AtomPair(atom2, atom3), distanceMatrix);
    }

    public QuickAndDirtyAngle(Angle angle) {
        super(angle.atomPair1,angle.atomPair2,angle.distanceMatrix);
    }
    /**
     * An approximated fast arc-cos function that uses the ArcCos class.
     * Note that the static method ArcCos.useFastArcCos() should be called
     * in order for the approximation to take affect. Standard Math.acos is the DEFAULT
     */
    public final double acos(double cos) {
        return ArcCos.acos(cos);
    }
}
