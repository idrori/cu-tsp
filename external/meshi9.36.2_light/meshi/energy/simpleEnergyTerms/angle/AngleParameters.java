/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.angle;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.geometry.*;
import meshi.util.filters.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.util.*;

import java.util.*;

public class AngleParameters implements Parameters, Comparable {
    public final double target, force, force2;
    public final AtomType type1, type2, type3;

    public AngleParameters() {
        this(AtomType.XXX, AtomType.XXX, AtomType.XXX, -9999.9, -9999.9);
    }

    public AngleParameters(AtomType type1, AtomType type2, AtomType type3) {
        this(type1, type2, type3, -9999.9, -9999.9);
    }

    public AngleParameters(String line) {
        this(new StringTokenizer(line));
    }

    private AngleParameters(StringTokenizer line) {
        this(AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                Utils.toDouble(line.nextToken()), Utils.toDouble(line.nextToken()));
    }

    public AngleParameters(AtomType type1, AtomType type2, AtomType type3, double targetAngle, double forceConstant) {
        this.type2 = type2;
        if (type1.compareTo(type3) < 0) {
            this.type1 = type1;
            this.type3 = type3;
        } else {
            this.type1 = type3;
            this.type3 = type1;
        }
        target = Angle.deg2rad(targetAngle);
        force = forceConstant;
        force2 = 2.0 * force;
    }

    public int compareTo(Object other) {
        if (!(other instanceof AngleParameters))
            throw new RuntimeException("Weird argument to " +
                    "AngleParameters.compairTo(Object other)");
        AngleParameters aep = (AngleParameters) other;
        if (type1.compareTo(aep.type1) > 0) return 1;
        if (type1.compareTo(aep.type1) < 0) return -1;
        if (type2.compareTo(aep.type2) > 0) return 1;
        if (type2.compareTo(aep.type2) < 0) return -1;
        if (type3.compareTo(aep.type3) > 0) return 1;
        if (type3.compareTo(aep.type3) < 0) return -1;
        return 0;
    }

    public String toString() {
        return "AngleParameters\n" +
                "\t type1          = " + type1 +
                "\t type2          = " + type2 +
                "\t type3          = " + type3 +
                "\t target   = " + Angle.rad2deg(target) + "\n" +
                "\t force = " + force;
    }

    public Parameters create(StringTokenizer line) {
        return new AngleParameters(line);
    }

    public Filter isA() {
        return new isA();
    }

    private class isA implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof AngleParameters);
        }
    }
}
