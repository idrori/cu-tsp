/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.PDB.PdbLine;
import meshi.util.*;
import meshi.molecularElements.atoms.*;

/**
 * --------------------------- Angle -----------------------------------
 */

public class Angle implements Updateable {
    public static final double EPSILON = 0.00000001;

    public static final double D5 = 0.0873; // ~5 degree  far too close to the numerically instable zero
    public static final double D175 =  3.0543; // ~175 degree  far too close to the numerically instable PI
    public static final double D10 = 0.1745; // ~10 degree  too close to the numerically instable zero
    public static final double D170 = 2.9671; // ~170 degree  too close to the numerically instable PI
    public static final double D160 = 2.7925;
    public static final double D150 =  2.6180;
    public static final double D1 =  0.0175;
    public static final double D179 = 3.1241;
    public static final double INSTABLE_LOW = D1;
    public static final double INSTABLE_HIGH = D179;
    public static final double TOO_HIGH = D160;
    public static final double TOO_LOW = D10;

    public final Atom atom1;
    public final Atom atom2;
    public final Atom atom3;
    protected String angleName = "";
    protected int angleCode = -1;
    protected int angleResNum = -9999;
    protected String angleResName = "";
    protected AtomPair atomPair1, atomPair2;
    public final Distance DISTANCE1, DISTANCE2;
    private double angle;
    private double dangleDx1, dangleDy1, dangleDz1;
    private double dangleDx2, dangleDy2, dangleDz2;
    private double dangleDx3, dangleDy3, dangleDz3;
    private double minusCosAngle;
    double  cosAngle;
    double  dCosAngleDx1;
    double  dCosAngleDx2;
    double  dCosAngleDx3;

    public double dCosAngleDx1() { return dCosAngleDx1; }
    public double dCosAngleDx2() { return dCosAngleDx2; }
    public double dCosAngleDx3() { return dCosAngleDx3; }

    public double dCosAngleDy1() { return dCosAngleDy1; }
    public double dCosAngleDy2() { return dCosAngleDy2; }
    public double dCosAngleDy3() { return dCosAngleDy3; }

    public double dCosAngleDz1() { return dCosAngleDz1; }
    public double dCosAngleDz2() { return dCosAngleDz2; }
    public double dCosAngleDz3() { return dCosAngleDz3; }

    double  dCosAngleDy1;
    double  dCosAngleDy2;
    double  dCosAngleDy3;
    double  dCosAngleDz1;
    double  dCosAngleDz2;
    double  dCosAngleDz3;

    private double sinAngle;
    protected DistanceMatrix distanceMatrix = null;
    private int numberOfUpdates = 0;

    /* Once upon a time there were BUG GENERATING constructors that did not use distance matrix.
       After they made enough damage they were removed */


    public Angle(Atom atom1, Atom atom2, Atom atom3, DistanceMatrix distanceMatrix) {
        this(new AtomPair(atom1, atom2), new AtomPair(atom2, atom3), distanceMatrix);
    }

    public Angle(AtomPair atomPair1, AtomPair atomPair2,
                 DistanceMatrix distanceMatrix) {
        this.atomPair1 = atomPair1;
        this.atomPair2 = atomPair2;
        this.distanceMatrix = distanceMatrix;
        Distance distance1, distance2;
        Distance distance1mirror, distance2mirror = null;

        if (atomPair1.atom1().number() > atomPair1.atom2().number()) {
            distance1 = distanceMatrix.distance(atomPair1.atom1(),
                    atomPair1.atom2());
            if (distance1 == null)
                throw new RuntimeException("Weird null distance between :\n"+
                                            atomPair1.atom1()+"\n"+atomPair1.atom1().core+"\n"+
                                            atomPair1.atom2()+"\n"+atomPair1.atom2().core+"\n"+
                                            "real distance = "+
                                            (new FreeDistance(atomPair1.atom1(),atomPair1.atom2())).distance());

            distance1mirror = new DistanceMirror(distance1);
        } else {
            distance1mirror = distanceMatrix.distance(atomPair1.atom2(),atomPair1.atom1());
            if (distance1mirror == null)
                throw new RuntimeException("Weird null distance between :\n"+
                                            atomPair1.atom1()+"\n"+atomPair1.atom1().core+"\n"+
                                            atomPair1.atom2()+"\n"+atomPair1.atom2().core+"\n"+
                                            "real distance = "+
                                            (new FreeDistance(atomPair1.atom1(),atomPair1.atom2())).distance());

            distance1 = new DistanceMirror(distance1mirror);

        }

        if (atomPair2.atom1().number() > atomPair2.atom2().number()) {
            distance2 = distanceMatrix.distance(atomPair2.atom1(),atomPair2.atom2());
            if (distance2 == null)
                throw new RuntimeException("Weird distance2 :\n"+
                        atomPair2.atom1()+"\n"+
                        atomPair2.atom1().core+"\n"+atomPair2.atom2()+"\n"+atomPair2.atom2().core);
                distance2mirror = new DistanceMirror(distance2);
        } else {
            distance2mirror = distanceMatrix.distance(atomPair2.atom2(),atomPair2.atom1());
            if (distance2mirror == null)
                throw new RuntimeException("Weird distance2 :\n"+
                        atomPair2.atom1()+"\n"+
                        atomPair2.atom1().core+"\n"+atomPair2.atom2()+"\n"+atomPair2.atom2().core);

            distance2 = new DistanceMirror(distance2mirror);
          }

        if (((atomPair1.atom1() == atomPair2.atom1()) && (atomPair1.atom2() == atomPair2.atom2())) ||
                ((atomPair1.atom1() == atomPair2.atom2()) && (atomPair1.atom2() == atomPair2.atom1())))
            throw new RuntimeException("Cannot create an angle from " +
                    atomPair1 +
                    " and " + atomPair2);

        if (atomPair1.atom1() == atomPair2.atom1()) {
            //  2-1
            //    1-2
            atom1 = atomPair1.atom2();
            atom2 = atomPair1.atom1();
            atom3 = atomPair2.atom2();
            DISTANCE1 = distance1mirror;
            DISTANCE2 = distance2mirror;
        } else if (atomPair1.atom1() == atomPair2.atom2()) {
            // 2-1
            //   2-1 
            atom1 = atomPair1.atom2();
            atom2 = atomPair1.atom1();
            atom3 = atomPair2.atom1();
            DISTANCE1 = distance1mirror;
            DISTANCE2 = distance2;
        } else if (atomPair1.atom2() == atomPair2.atom2()) {
            // 1-2
            //   2-1
            atom1 = atomPair1.atom1();
            atom2 = atomPair1.atom2();
            atom3 = atomPair2.atom1();
            DISTANCE1 = distance1;
            DISTANCE2 = distance2;
        } else if (atomPair1.atom2() == atomPair2.atom1()) {
            // 1-2
            //   1-2
            atom1 = atomPair1.atom1();
            atom2 = atomPair1.atom2();
            atom3 = atomPair2.atom2();
            DISTANCE1 = distance1;
            DISTANCE2 = distance2mirror;
        } else throw new RuntimeException("Cannot create an angle from " +
                atomPair1 +
                " and " + atomPair2);
        assignName();
        try {
            update();
        } catch (AngleException ae) {
        } // This exception is thrown when the angle is too close to zero or PI. By ignoring it here we allow the calling method to decide what to do with this angle.

    }

    public Atom atom1() {
        return atom1;
    }

    public Atom atom2() {
        return atom2;
    }

    public Atom atom3() {
        return atom3;
    }

    public AtomPair atomPair1() {
        return atomPair1;
    }

    public AtomPair atomPair2() {
        return atomPair2;
    }

    public AtomPair sharedAtomPair(Angle other) {
        if (atomPair1.equals(other.atomPair1)) return atomPair1;
        if (atomPair1.equals(other.atomPair2)) return atomPair1;
        if (atomPair2.equals(other.atomPair1)) return atomPair2;
        if (atomPair2.equals(other.atomPair2)) return atomPair2;
        return null;
    }

    public boolean proper(Angle other) {
        AtomPair shared = sharedAtomPair(other);
        return shared != null && atom2 != other.atom2;
    }

    public AtomPair hinge(Angle other) {
        AtomPair out = sharedAtomPair(other);
        if (out == null) return null;
        if ((atom2 == other.atom1()) |
                (atom2 == other.atom3())) return out;
        return null;
    }

    public String toString() {
        return "Angle: " + angleName +
                "\n\tatom1   :" + atom1() + " " + atom1().getClass() +
                "\n\tatom2   :" + atom2() + " " + atom2().getClass() +
                "\n\tatom3   :" + atom3() + " " + atom3().getClass() +
                "\n\tvalue   :" + rad2deg(angle());
    }

    public double angle() {
        return angle;
    }

    public double dangleDx1() {
        return dangleDx1;
    }

    public double dangleDy1() {
        return dangleDy1;
    }

    public double dangleDz1() {
        return dangleDz1;
    }

    public double dangleDx2() {
        return dangleDx2;
    }

    public double dangleDy2() {
        return dangleDy2;
    }

    public double dangleDz2() {
        return dangleDz2;
    }

    public double dangleDx3() {
        return dangleDx3;
    }

    public double dangleDy3() {
        return dangleDy3;
    }

    public double dangleDz3() {
        return dangleDz3;
    }

    public double cosAngle() {
        return cosAngle;
    } // Bug fixed by Chen 6.5.2019

    public double sinAngle() {
        return sinAngle;
    }

    public void update(int numberOfUpdates) throws UpdateableException {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            update();
            this.numberOfUpdates++;
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with Angle.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }

    public void update() {
        double dd1Dx = DISTANCE1.dDistanceDx();
        double dd1Dy = DISTANCE1.dDistanceDy();
        double dd1Dz = DISTANCE1.dDistanceDz();
        double invD1 = DISTANCE1.invDistance();
        double dd2Dx = -1 * DISTANCE2.dDistanceDx();
        double dd2Dy = -1 * DISTANCE2.dDistanceDy();
        double dd2Dz = -1 * DISTANCE2.dDistanceDz();
        double invD2 = DISTANCE2.invDistance();

        minusCosAngle = (dd1Dx * dd2Dx + dd1Dy * dd2Dy + dd1Dz * dd2Dz);
        cosAngle = -minusCosAngle;

        double d_dd1Dx_dx =  invD1*(1-dd1Dx*dd1Dx);
        double d_dd1Dx_dy = -invD1*dd1Dy*dd1Dx;
        double d_dd1Dx_dz = -invD1*dd1Dz*dd1Dx;

        double d_dd1Dy_dx = -invD1*dd1Dy*dd1Dx;
        double d_dd1Dy_dy =  invD1*(1 - dd1Dy*dd1Dy);
        double d_dd1Dy_dz = -invD1*dd1Dz*dd1Dy;

        double d_dd1Dz_dx = -invD1*dd1Dx*dd1Dz;
        double d_dd1Dz_dy = -invD1*dd1Dy*dd1Dz;
        double d_dd1Dz_dz =  invD1*(1 - dd1Dz*dd1Dz);

        dCosAngleDx1 = -dd2Dx * d_dd1Dx_dx - dd2Dy * d_dd1Dx_dy - dd2Dz * d_dd1Dx_dz ;
        dCosAngleDy1 = -dd2Dy * d_dd1Dy_dy - dd2Dx * d_dd1Dy_dx - dd2Dz * d_dd1Dy_dz ;
        dCosAngleDz1 = -dd2Dz * d_dd1Dz_dz - dd2Dy * d_dd1Dz_dy - dd2Dx * d_dd1Dz_dx ;

        double d_dd2Dx_dx = -invD2*(1 - dd2Dx*dd2Dx);
        double d_dd2Dx_dy =  invD2*dd2Dy*dd2Dx;
        double d_dd2Dx_dz =  invD2*dd2Dz*dd2Dx;

        double d_dd2Dy_dx =  invD2*dd2Dy*dd2Dx;
        double d_dd2Dy_dy = -invD2*(1 - dd2Dy*dd2Dy);
        double d_dd2Dy_dz =  invD2*dd2Dz*dd2Dy;

        double d_dd2Dz_dx =  invD2*dd2Dx*dd2Dz;
        double d_dd2Dz_dy =  invD2*dd2Dy*dd2Dz;
        double d_dd2Dz_dz = -invD2*(1 - dd2Dz*dd2Dz);

        dCosAngleDx3 = -dd1Dx * d_dd2Dx_dx - dd1Dy * d_dd2Dx_dy - dd1Dz * d_dd2Dx_dz ;
        dCosAngleDy3 = -dd1Dy * d_dd2Dy_dy - dd1Dx * d_dd2Dy_dx - dd1Dz * d_dd2Dy_dz ;
        dCosAngleDz3 = -dd1Dz * d_dd2Dz_dz - dd1Dy * d_dd2Dz_dy - dd1Dx * d_dd2Dz_dx ;

        dCosAngleDx2 = -dCosAngleDx1 - dCosAngleDx3;
        dCosAngleDy2 = -dCosAngleDy1 - dCosAngleDy3;
        dCosAngleDz2 = -dCosAngleDz1 - dCosAngleDz3;


        //The derivatives and torsion angles become unstable when the the angle gets too close to zero or PI.
        // One may of course ignore this exception if derivation is irelevant or handled some other way.
        if ((minusCosAngle >= 1.0 - EPSILON) || (minusCosAngle <= -1.0 + EPSILON)) {
            if (minusCosAngle >= 1.0 - EPSILON) angle = Math.PI;
            else angle = 0;
//            throw new AngleException(this, minusCosAngle);
        }
        else {
        angle = acos(-1 * minusCosAngle);
        double sinSign;
        if ((angle >= 0) & (angle <= Math.PI))
            sinSign = 1;
        else
            sinSign = -1;
        sinAngle = sinSign*Math.sqrt(1 - cosAngle*cosAngle);
        //sinAngle = Math.sin(angle);
        double factor = invD1 / sinAngle;
        dangleDx1 = factor * (dd2Dx - minusCosAngle * dd1Dx);
        dangleDy1 = factor * (dd2Dy - minusCosAngle * dd1Dy);
        dangleDz1 = factor * (dd2Dz - minusCosAngle * dd1Dz);
        factor = invD2 / sinAngle;
        dangleDx3 = factor * (minusCosAngle * dd2Dx - dd1Dx);
        dangleDy3 = factor * (minusCosAngle * dd2Dy - dd1Dy);
        dangleDz3 = factor * (minusCosAngle * dd2Dz - dd1Dz);
        dangleDx2 = -1 * (dangleDx1 + dangleDx3);
        dangleDy2 = -1 * (dangleDy1 + dangleDy3);
        dangleDz2 = -1 * (dangleDz1 + dangleDz3);
        }
    }

    /**
     * Calculate the arc-cosine function.
     * This method simply calls Math.acos, extending classes however may
     * use approximations in order to save time.
     *
     * @param cos value
     * @return arccos value
     */
    public double acos(double cos) {
        return Math.acos(cos);
    }

    public boolean frozen() {
        return (atom1.frozen() & atom2.frozen() & atom3.frozen());
    }

    public static double deg2rad(double ang) {
        return ang * Math.PI / 180;
    }

    public static double rad2deg(double ang) {
        return ang * 180 / Math.PI;
    }


    public String getAngleName() {
        return angleName;
    }

    public int getAngleCode() {
        return angleCode;
    }

    public int getAngleResNum() {
        return angleResNum;
    }

    public String getAngleResName() {
        return angleResName;
    }

    /**
     * Assigns a meaningful name to an angle.
     * This method assigns (if possible):
     * - name to the angle object.
     * - number corresponding to this name.
     * - the number in the chain of the center atom residue to which this amgle belong
     * - the name of the center atom residue to which this angle belong
     * <p/>
     * The names to numbers conversion is:
     * Three consecutive Ca's - 0
     */
    protected void assignName() {
        int resNum = atom2.residueNumber();
        if ((atom1.name().compareTo("CA") == 0) &&                 // CA3
                (atom2.name().compareTo("CA") == 0) &&
                (atom3.name().compareTo("CA") == 0) &&
                (((atom3.residueNumber() == (resNum + 1)) && (atom1.residueNumber() == (resNum - 1))) ||
                        ((atom1.residueNumber() == (resNum + 1)) && (atom3.residueNumber() == (resNum - 1))))) {
            angleName = "CA3";
            angleCode = 0;
            angleResNum = resNum;
            angleResName = atom2.residueName();
        }
    }

    public static class AngleException extends RuntimeException {
        private static String errorMessage = "Weird angle:\n%s\n" + "with minusCosAngle = %16.14f while it should be between %16.15f and %16.15f.\n" +
                "The derivatives and torsion angles become unstable when the the angle gets too close to zero or PI.\n" +
                "One may of course ignore this exception if derivation is irelevant or handled some other way.";

        public AngleException(double cosAngle) {
            super(String.format(errorMessage, "Unknown angle", cosAngle, -1 + EPSILON, 1 - EPSILON));
        }

        public AngleException(Angle angle, double cosAngle) {
            super(String.format(errorMessage, angle.toString(), cosAngle, -1 + EPSILON, 1 - EPSILON));
        }
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }

    public boolean dangerous() {
        if ((angle > INSTABLE_HIGH) || (angle < INSTABLE_LOW)) return true;
        return false;
    }
    public boolean disturbing() {
        if ((angle > TOO_HIGH) || (angle < TOO_LOW)) return true;
        return false;
    }

    public static void main(String[] args) throws UpdateableException{
        double epsilon = 0.0000000001;

        PdbLine pdbLine1 = new PdbLine("ATOM      1  N   MET     1       0.000   0.000   0.000  1.00  0.00");
        PdbLine pdbLine2 = new PdbLine("ATOM      0  CA  MET     1       0.000   1.330   0.000  1.00  0.00");
        PdbLine pdbLine3 = new PdbLine("ATOM      3  C   MET     1       0.000   1.880  -1.417  1.00  0.00");

        MolecularSystem ms = new MolecularSystem();
        Atom atom1 = new Atom(pdbLine1, ms);
        Atom atom2 = new Atom(pdbLine2, ms);
        Atom atom3 = new Atom(pdbLine3, ms);
        ms.terminator().reset();
        DistanceMatrix distanceMatrix = new DistanceMatrix(ms);
        Angle angle = new Angle(atom1, atom2, atom3, distanceMatrix);
        System.out.println(angle.atom1);
        System.out.println(angle.atom2);
        System.out.println(angle.atom3);
        System.out.println(angle.DISTANCE1.atom1);
        System.out.println(angle.DISTANCE1.atom2);
        System.out.println(angle.DISTANCE2.atom1);
        System.out.println(angle.DISTANCE2.atom2);

        distanceMatrix.update();
        angle.update();
        double c1 = angle.cosAngle();
        double a = angle.dCosAngleDz2();
        atom2.addToZ(epsilon);
        distanceMatrix.update();
        angle.update();
        double c2 = angle.cosAngle();
        System.out.println(angle.angle()+" "+angle.cosAngle+" "+a+" "+(c2-c1)/epsilon);
    }

}
