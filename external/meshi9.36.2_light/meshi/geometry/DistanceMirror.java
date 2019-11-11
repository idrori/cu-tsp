/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

import meshi.util.*;

public class DistanceMirror extends Distance {
    public final Distance source;

    public DistanceMirror(Distance source) {
        super(source.atom2, source.atom1, -999999.999, -999999.999, -999999.999, -999999.999);
        this.source = source;
        if (atom2Number < atom1Number) throw new RuntimeException("weird distance mirror \n" +
                "atom2Number = " + atom2Number + "  tom1Number = " + atom1Number + "\n" +
                atom1() + "\n" + atom2());
        mode = DistanceMode.MIRROR;
    }

    public final double distance() {
        return source.distance();
    }

    public final double invDistance() {
        return source.invDistance();
    }

    public final double dDistanceDx() {
        return 0 - source.dDistanceDx();
    }

    public final double dDistanceDy() {
        return 0 - source.dDistanceDy();
    }

    public final double dDistanceDz() {
        return 0 - source.dDistanceDz();
    }

    public final double dx() {
        return 0 - source.dx();
    }

    public final double dy() {
        return 0 - source.dy();
    }

    public final double dz() {
        return 0 - source.dz();
    }

    public final boolean dead() {
        return source.dead();
    }

}

