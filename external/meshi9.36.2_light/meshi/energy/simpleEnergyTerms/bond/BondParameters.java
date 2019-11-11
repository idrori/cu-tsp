/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.bond;

import meshi.energy.*;
import meshi.parameters.*;
import meshi.util.*;
import meshi.util.filters.*;

import java.util.*;

/**
 * The energy parameters associated with a bond between two getAtoms.
 * These parameters depend isOn the atom types.
 */
public class BondParameters implements Parameters, Comparable {
    /**
     * The target distance for this type of bond.
     */
    public final double target;

    /**
     * The force constant for this type of bond.
     */
    public final double force;

    /**
     * The square of the force constant.
     * Not a real parameter but saves time.
     */
    public final double force2;

    /**
     * Atom type of the first atom.
     * The convention is that it is the smaller atom type in the pair.
     */
    public final AtomType type1;

    /**
     * Atom type of the second atom.
     * The convention is that it is the larger atom type in the pair.
     */
    public final AtomType type2;

    /**
     * Auxiliary constructor.
     */
    public BondParameters() {
        this(AtomType.XXX, AtomType.XXX, -999, 999);
    }

    /**
     * Constructs BondParmeters from a string  (typically a line in the bondParameters file).
     */
    public BondParameters(String line) {
        this(new StringTokenizer(line));
    }

    /**
     * A helper constructor for ( BondParameters(String line) ).
     */
    public BondParameters(StringTokenizer line) {
        this(AtomType.type(line.nextToken()), // type1
                AtomType.type(line.nextToken()), // type2
                Utils.toDouble(line.nextToken()), // targetDistance
                Utils.toDouble(line.nextToken())); // forceConstant
    }

    /**
     * A helper constructor for ( BondParameters(String line) ).
     */
    public BondParameters(AtomType type1, AtomType type2,
                          double targetDistance, double forceConstant) {
        if (type1.compareTo(type2) > 0) {
            this.type1 = type1;
            this.type2 = type2;
        } else {
            this.type1 = type2;
            this.type2 = type1;
        }
        target = targetDistance;
        force = forceConstant;
        force2 = 2.0 * forceConstant;
    }

    /**
     * For search keys in parameterList
     */
    public BondParameters(AtomType type1, AtomType type2) {
        this(type1, type2, -1, -1);
    }


    /**
     * Defines order within objects and thus allows sorting for more efficient searches.
     */
    public int compareTo(Object other) {
        if (!(other instanceof BondParameters))
            throw new RuntimeException("Weird argument to " +
                    "BondParameters.compairTo(Object other)");
        BondParameters bp = (BondParameters) other;
        if (type1.compareTo(bp.type1) > 0) return 1;
        if (type1.compareTo(bp.type1) < 0) return -1;
        if (type2.compareTo(bp.type2) > 0) return 1;
        if (type2.compareTo(bp.type2) < 0) return -1;
        return 0;
    }

    public String toString() {
        return "BondParameters\n" +
                "\t type1  = " + type1 + "\n" +
                "\t type2  = " + type2 + "\n" +
                "\t target = " + target + "\n" +
                "\t force  = " + force;
    }

    private class isA implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof BondParameters);
        }
    }
}
