/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.*;
import meshi.util.formats.*;
import meshi.util.*;
import meshi.parameters.*;
import meshi.util.Attributable;
import meshi.util.MeshiAttribute;

/**
 * The distance between two {@link meshi.molecularElements.atoms}<br>.
 * <p/>
 * Almost any measurable feature of a molecule is related to distances
 * between pairs of {@link meshi.molecularElements.atoms}.
 * Thus, The calculation of distances (and their inverse and derivatives)
 * are typically a computational bottleneck in computational structural
 * biology applications. In all applications that we are aware of (Please,
 * enlighten us if you know better) distance calculation is done as part
 * of the procedures that use it (say, as part of the van-der-Waals energy
 * calculation). As a result the distance between two getAtoms may be calculated
 * more then once. For example the distance between two getAtoms may be
 * calculated both during angle and torsion angle energies calculations.
 * In Meshi We tried to consentrate all distance related issues in a few
 * classes: this one, its subclasses and the closely connected
 * class  {@link DistanceMatrix DistanceMatrix}.
 * <p/>
 * <b> Possible pitfalls </b><br>
 * <ol>
 * <li>The object variables dx, dy, dz, distance, distance2, invDistance,
 * dDistanceDx, dDistanceDy and dDistanceDz should have been private and
 * accessed through "get methods". Energy functions though, use the
 * values of these variables intensively and a considerable computational
 * gain is achieved granting them public accessability (protected for
 * dx, dy and dz) and removing the function call overhead. Thus,
 * changing the values of these variables by other classes is
 * possible (from the compiler point of view) but not very much
 * recommended.
 * <li> The distance object <b>"is not aware" </b> of changes in the
 * getAtoms coordinates. The values stored in this object are correct
 * only after the {@link meshi.geometry.Distance#update update}
 * method is <u>explicitly</u> called.
 */
public class Distance implements Attributable {
    //---------------------------- fields ------------------------------
    public static final double INFINITE_DISTANCE = Double.MAX_VALUE;
    /**
     * This object represent the distance between atom1 and atom2.
     */
    public final AtomCore atom1, atom2;
    public final int atom1Number;
    protected int atom2Number;
    protected DistanceMode mode;
    public final boolean heavy;

    /**
     * The inverse of the distance between atom1 & atom2.         <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door
     * for bugs.
     */
    public double invDistance;


    /**
     * The distance between atom1 & atom2  <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door
     * for bugs.
     */
    public double distance;

    /**
     * atom1.x - atom2.x
     * Not that the protected accessability of of this variable
     * improves computational efficiency but opens a wide door
     * for bugs.
     */
    public double dx;
    /**
     * atom1.y - atom2.y
     * Not that the protected accessability of of this variable
     * improves computational efficiency but opens a wide door
     * for bugs.
     */
    public double dy;
    /**
     * atom1.z - atom2.z
     * Not that the protected accessability of of this variable
     * improves computational efficiency but opens a wide door
     * for bugs.
     */
    public double dz;

    public final AtomType largeType, smallType;
    private AttributesRack attributes = new AttributesRack();

    public final void addAttribute(MeshiAttribute attribute) {
        attributes.addAttribute(attribute);
    }

    public final MeshiAttribute getAttribute(int key) {
        return attributes.getAttribute(key);
    }

    //-------------------------- Constructors --------------------------

    public Distance(AtomCore atom1, AtomCore atom2, double dx, double dy, double dz, double distance) {
        this.atom1 = atom1;
        this.atom2 = atom2;
        atom1Number = atom1.number;
        atom2Number = atom2.number;
        this.dx = dx;
        this.dy = dy;
        this.dz = dz;
        this.distance = distance;
        invDistance = 1 / distance;
        mode = null;
        if (atom1.type().compareTo(atom2.type()) > 0) {
            largeType = atom1.type();
            smallType = atom2.type();
        } else {
            largeType = atom2.type();
            smallType = atom1.type();
        }
        if (atom1.type().isHydrogen() || atom2.type().isHydrogen()) heavy = false;
        else heavy = true;
    }

    public void update(double rMax2) {
        if ((atom1 == null) || (atom2 == null))
            throw new RuntimeException("Try to update weird distance "+this);
        dx = atom1.x() - atom2.x();
        dy = atom1.y() - atom2.y();
        dz = atom1.z() - atom2.z();
        double distance2 = dx * dx + dy * dy + dz * dz;
        if ((distance2 < rMax2) || this.mode.bonded) {
            distance = Math.sqrt(distance2);
            if (distance < 0.000001) throw new RuntimeException("Distance too close to zero " + this);
            invDistance = 1 / distance;
        }
        else{
            mode = DistanceMode.INFINITE;
            distance = Distance.INFINITE_DISTANCE;
        }
    }


    //------------------------------- methods ------------------------------

    /**
     * Get the distance between atom1 and atom2.
     */
    public double distance() {
        return distance;
    }

    /**
     * Get the inverse of the distance between atom1 and atom2.
     */
    public double invDistance() {
        return invDistance;
    }


    /**
     * The derivative of the distance between atom1 & atom2 by
     * the X coordinate of atom1. The derivative by the X coordinate of
     * atom2 is this value multiplied by -1.
     */
    public double dDistanceDx() {
        return dx * invDistance;
    }

    /**
     * The derivative of the distance between atom1 & atom2 by
     * the Y coordinate of atom1. The derivative by the Y coordinate of
     * atom2 is this value multiplied by -1.
     */
    public double dDistanceDy() {
        return dy * invDistance;
    }

    /**
     * The derivative of the distance between atom1 & atom2 by
     * the Z coordinate of atom1. The derivative by the Z coordinate of
     * atom2 is this value multiplied by -1.
     */
    public double dDistanceDz() {
        return dz * invDistance;
    }

    /**
     * Get atom1.
     */
    public Atom atom1() {
        return atom1.atom;
    }

    /**
     * Get atom2.
     */
    public Atom atom2() {
        return atom2.atom;
    }

    public int atom2Number() {
        return atom2Number;
    }

    public double dx() {
        return dx;
    }

    public double dy() {
        return dy;
    }

    public double dz() {
        return dz;
    }

    public String toString() {
        Fdouble dformat = Fdouble.STANDARD;
        return "" + this.getClass() + " " + mode + " " + atom1.number + " " + atom1.atom.residueNumber() + " " + atom2.number + " " + atom2.atom.residueNumber() + " " +
                dformat.f(distance()) + " " +
                dformat.f(dDistanceDx()) + " " +
                dformat.f(dDistanceDy()) + " " +
                dformat.f(dDistanceDz());
    }

    public double summaValue = 0;      //Todo protect this field

    public final DistanceMode mode() {
        return mode;
    }

    public void setMode(DistanceMode mode) {
        this.mode = mode;
    }

    public boolean dead() {
        return mode.dead;
    }


}

