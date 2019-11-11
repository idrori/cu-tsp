/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.plane;

import meshi.energy.Parameters;
import meshi.parameters.AtomType;
import meshi.util.Utils;
import meshi.util.filters.Filter;

import java.util.StringTokenizer;

public class PlaneParameters implements Parameters, Comparable {
    public enum Isomer {TRANS, CIS, CIS_TRANS;}
    /*public final static int TRANS = 1;
    public final static int CIS = 2;
    public final static int CIS_TRANS = 3;
*/
    public final double force, force2;
    public Isomer isomer;
    public final AtomType type1, type2, type3, type4;
    static String temp;

    public PlaneParameters() {
        this(AtomType.XXX, AtomType.XXX, AtomType.XXX, AtomType.XXX, -9999.9, Isomer.CIS_TRANS);
    }

    public PlaneParameters(AtomType type1, AtomType type2, AtomType type3, AtomType type4) {
        this(type1, type2, type3, type4, -9999.9, Isomer.CIS_TRANS);
    }

    public PlaneParameters(String line) {
        this(new StringTokenizer(line));
    }

    private PlaneParameters(StringTokenizer line) {
        this(AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                AtomType.type(line.nextToken()),
                Utils.toDouble(line.nextToken()),
                getType(line.nextToken()));
    }

    private static Isomer getType(String token) {
        if (token.equals("TRANS")) return Isomer.TRANS;
        if (token.equals("CIS")) return Isomer.CIS;
        return Isomer.CIS_TRANS;
    }
    public PlaneParameters(AtomType type1, AtomType type2,
                           AtomType type3, AtomType type4,
                           double force, Isomer isomer) {
        this.type1 = type1;
        this.type2 = type2;
        this.type3 = type3;
        this.type4 = type4;
        this.force = force;
        force2 = 2 * force;
        this.isomer = isomer;
    }

    public int compareTo(Object obj) {
        PlaneParameters other = (PlaneParameters) obj;
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
        return "PlaneParameters\n" +
                "\t type1   = " + type1 +
                "\t type2   = " + type2 +
                "\t type3   = " + type3 +
                "\t type4   = " + type4 +
                "\t isomer   = " + isomer + "\n" +
                "\t force = " + force;
    }

    public Parameters create(StringTokenizer line) {
        return new PlaneParameters(line);
    }

    public Filter isA() {
        return new isA();
    }

    private class isA implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof PlaneParameters);
        }
    }
}
