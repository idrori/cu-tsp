/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

/**
 * Allows the use of approximated arc-cos function.
 * The use of the standard Math.acos function is the default.
 */

public class QuickAndDirtyTorsion extends Torsion {
    //private final ArcCos ARC_COS = new ArcCos();

    public QuickAndDirtyTorsion(Angle angle1, Angle angle2, DistanceMatrix distanceMatrix) {
        super(new QuickAndDirtyAngle(angle1), new QuickAndDirtyAngle(angle2), distanceMatrix);
    }

    /**
     * An approximated fast arc-cos function that uses the ArcCos class.
     * Note that the static method ArcCos.useFastArcCos() should be called
     * in order for the approximation to take affect. Standard Math.acos is the DEFAULT
     */
    public double acos(double cos) {
        return ArcCos.acos(cos);
    }
}
