/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.*;

public class BondedDistance extends Distance {
    public BondedDistance(Atom atom1, Atom atom2, double dx, double dy, double dz, double distance) {
        super(atom1.core, atom2.core, dx, dy, dz, distance);
        mode = DistanceMode.BONDED;
    }

    protected void update() {
        if ((atom1 == null) || (atom2 == null))
            throw new RuntimeException("Try to update weird distance");
        dx = atom1.x() - atom2.x();
        dy = atom1.y() - atom2.y();
        dz = atom1.z() - atom2.z();
        distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        invDistance = 1 / distance;
    }
}

