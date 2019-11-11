/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

/**
 * A list of quick and dirty torsion angles. These angle may use approximate
 * torsion calculations if ArcCos.useFastArcCos() is called.
 */
public class QuickAndDirtyTorsionList extends TorsionList {
    public QuickAndDirtyTorsionList(AngleList angles, DistanceMatrix distanceMatrix) {
        super(angles, distanceMatrix);
    }

    public Torsion getTorsion(Angle angle1, Angle angle2, DistanceMatrix distanceMatrix) {
        return new QuickAndDirtyTorsion(angle1, angle2, distanceMatrix);
    }
}
