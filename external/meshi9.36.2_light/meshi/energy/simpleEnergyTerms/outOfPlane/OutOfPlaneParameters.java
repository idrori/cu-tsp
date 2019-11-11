/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.outOfPlane;

import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.geometry.*;
import meshi.util.filters.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.*;

import java.util.*;

public class OutOfPlaneParameters implements Parameters, Comparable {
    public final double target, force, force2;
    public final AtomType type1, type2, type3, type4;

    public OutOfPlaneParameters() {
        this(AtomType.XXX, AtomType.XXX, AtomType.XXX, AtomType.XXX, -9999.9, -9999.9);
    }

    public OutOfPlaneParameters(AtomType type1, AtomType type2, AtomType type3, AtomType type4) {
        this(type1, type2, type3, type4, -9999.9, -9999.9);
    }

    public OutOfPlaneParameters(String line) {
        this(new StringTokenizer(line));
    }

    private OutOfPlaneParameters(StringTokenizer line) {
        this(AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                Utils.toDouble(line.nextToken()), Utils.toDouble(line.nextToken()));
    }

    public OutOfPlaneParameters(AtomType type1, AtomType type2,
                                AtomType type3, AtomType type4,
                                double target, double force) {
        this.type1 = type1;
        this.type2 = type2;
        this.type3 = type3;
        this.type4 = type4;
        this.force = force;
        force2 = 2 * force;
        this.target = Angle.deg2rad(target);
    }

    public int compareTo(Object obj) {
        OutOfPlaneParameters other = (OutOfPlaneParameters) obj;
        if (type1.compareTo(other.type1) > 0) return 1;
        if (type1.compareTo(other.type1) < 0) return -1;
        if (type2.compareTo(other.type2) > 0) return 1;
        if (type2.compareTo(other.type2) < 0) return -1;
        if (type3.compareTo(other.type3) > 0) return 1;
        if (type3.compareTo(other.type3) < 0) return -1;
        if (type4.compareTo(other.type4) > 0) return 1;
        if (type4.compareTo(other.type4) < 0) return -1;
        return 0;
    }

    public String toString() {
        return "OutOfPlaneParameters\n" +
                "\t type1          = " + type1 +
                "\t type2          = " + type2 +
                "\t type3          = " + type3 +
                "\t type4          = " + type4 +
                "\t target   = " + Angle.rad2deg(target) + "\n" +
                "\t force = " + force;
    }

    public Parameters create(StringTokenizer line) {
        return new OutOfPlaneParameters(line);
    }

    public Filter isA() {
        return new isA();
    }

    private class isA implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof OutOfPlaneParameters);
        }
    }
}
