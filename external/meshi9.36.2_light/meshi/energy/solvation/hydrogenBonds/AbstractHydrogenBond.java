/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation.hydrogenBonds;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.*;


/**
 * A super class with the most general treatment of hydrogen bond.
 * A specific implementation should put content to these abstract methods.
 * The constructor gets a distance matrix object, and an AtomList object of
 * all the getAtoms participating in the bond. The only condition isOn the list is
 * that the donor and acceptor must be its first two elements (their relative
 * order is not important).
 * <p/>
 * IPORTANT IMPORTANT IMPORTANT
 * NOTE: The distance matrix must be in an updated state before any of this class's
 * methods are called. If it is not up to date, the results would be meaningless!!!
 * This class does not update the distance matrix isOn its own at any time!!
 */
public abstract class AbstractHydrogenBond {

    protected AtomList atomList = null;
    protected DistanceMatrix dm = null;
    protected double hbVal = 0.0;

    /**
     * derivatives[i][0] is the derivative of the HB value relative to the X coordinate of atom i.
     * derivatives[i][1] is the derivative of the HB value relative to the Y coordinate of atom i.
     * derivatives[i][2] is the derivative of the HB value relative to the Z coordinate of atom i.
     */
    protected double[][] derivatives = null;

    public AbstractHydrogenBond() {
    }

    public AbstractHydrogenBond(DistanceMatrix dm, AtomList atomList) {
        this.dm = dm;
        this.atomList = atomList;
        derivatives = new double[atomList.size()][3];
    }


    /**
     * This update method should be used if a change to the getAtoms coordinates has occured.
     * Only after this update, the getter methods can produce correct values.
     * This class assume the distance matrix is updated.
     */
    public abstract void updateHBvalueAndDerivatives();


    /**
     * FOLLOWING the activation of "update" this method can be activated,
     * where it will apply the forces to the getAtoms participating in the HB.
     * The forces are: factor*(minus derivative). This class is useful in energy terms that
     * make use of the hydrogen bond.
     */
    public void applyForcesToAtoms(double factor) {
        Atom atom;

        for (int c = 0; c < atomList.size(); c++) {
            atom = atomList.atomAt(c);
//            if (derivatives[c][0] + derivatives[c][1] + derivatives[c][2]> 50)
//                System.out.println("Hot hot hot "+derivatives[c][0] + derivatives[c][1] + derivatives[c][2]+
//                        "\n"+atom);
            if (!atom.frozen()) {
                atomList.atomAt(c).addToFx(-factor * derivatives[c][0]);
                atomList.atomAt(c).addToFy(-factor * derivatives[c][1]);
                atomList.atomAt(c).addToFz(-factor * derivatives[c][2]);
            }
        }
    }



    public double hbVal() {
        return hbVal;
    }

    public double derivatives(int atomInd, int xyz) {
        return derivatives[atomInd][xyz];
    }

    public final double[][] derivatives() {
        return derivatives;
    }

    public final AtomList atomList() {
        return atomList;
    }

    public final Atom getFirstPolar() {
        return atomList.atomAt(0);
    }

    public final Atom getSecondPolar() {
        return atomList.atomAt(1);
    }

    public String toString() {
        return "HB: " + atomList.atomAt(0).residue() + " " + atomList.atomAt(0).name() + " --- " +
                atomList.atomAt(1).residue() + " " + atomList.atomAt(1).name();
	}
}


  
