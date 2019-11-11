/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.*;

public class BondedFrozenDistance extends FrozenDistance {
    public BondedFrozenDistance(Atom atom1, Atom atom2, double dx, double dy, double dz, double distance) {
        super(atom1.core, atom2.core, dx, dy, dz, distance);
        mode = DistanceMode.BONDED_FROZEN;
    }
}

