/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.filters.*;
import meshi.util.*;
import meshi.parameters.*;

/**
 * The distance between two {@link meshi.molecularElements.Atom Atoms}<br>.
 * <p/>
 * Almost any measurable feature of a molecule is related to distances
 * between pairs of {@link meshi.molecularElements.Atom Atoms}.
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
 * only after the {@link Distance#update update}
 * method is <u>explicitly</u> called.
 */
public class FreeDistance extends Distance {
    private static double d2;
    //-------------------------- Constructors --------------------------

    public FreeDistance(Atom atom1, Atom atom2) {
        super(atom1.core, atom2.core, getDx(atom1, atom2), getDy(atom1, atom2), getDz(atom1, atom2), getDistance(atom1, atom2));
        mode = DistanceMode.FREE;
        update();
    }
    //------------------------------- methods ------------------------------


    public static double getDx(Atom atom1, Atom atom2) {
        return atom1.x() - atom2.x();
    }

    public static double getDy(Atom atom1, Atom atom2) {
        return atom1.y() - atom2.y();
    }

    public static double getDz(Atom atom1, Atom atom2) {
        return atom1.z() - atom2.z();
    }

    public static double getDistance(Atom atom1, Atom atom2) {
        double dx = atom1.x() - atom2.x();
        double dy = atom1.y() - atom2.y();
        double dz = atom1.z() - atom2.z();
        double d2 = dx * dx + dy * dy + dz * dz;
        return Math.sqrt(d2);
    }

    public void update() {
        dx = atom1.x() - atom2.x();
        dy = atom1.y() - atom2.y();
        dz = atom1.z() - atom2.z();
        d2 = dx * dx + dy * dy + dz * dz;
        distance = Math.sqrt(d2);
//        if (distance < 0.000001) throw new RuntimeException("Distance too close to zero "+this);

        invDistance = 1 / distance;
    }
}

