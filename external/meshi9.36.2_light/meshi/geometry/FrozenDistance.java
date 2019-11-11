/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.*;

public class FrozenDistance extends Distance {
    public FrozenDistance(AtomCore atom1, AtomCore atom2, double dx, double dy, double dz, double distance) {
        super(atom1, atom2, dx, dy, dz, distance);
        mode = DistanceMode.FROZEN;
    }

}

