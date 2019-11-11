/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.inflate.inflate;

import meshi.energy.EnergyElement;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.MeshiProgram;
//-------------------------------------------------------------------------------------------------------

//                                     InflateElement
//-------------------------------------------------------------------------------------------------------
public class InflateEnergyElement extends EnergyElement {
    protected Atom atom1, atom2;
    //    protected AtomPair atomPair;
    protected Atom atom1Copy, atom2Copy;
    protected double dEdD;
    protected double dEdX;
    protected double dEdY;
    protected double dEdZ;
    protected double energy;
    protected double weight;
    protected FreeDistance distance;
    protected double target;
    public static final double ALPHA = 0.1;
    private double delta;

    public InflateEnergyElement(FreeDistance distance, double weight) {
        this(distance, distance.distance(), 10.0, weight);
    }

    public InflateEnergyElement(FreeDistance distance, double targetDistance, double range, double weight) {
        atom1 = distance.atom1();
        atom2 = distance.atom2();
        setAtoms();
        updateFrozen();
        this.distance = distance;
        this.weight = weight;
        delta = range * (MeshiProgram.randomNumberGenerator().nextDouble() - 0.5);
        distance.update();
        target = targetDistance + delta;
    }

    public String toString() {
        return "RamachInflateElementOf " + distance + " delta = " + delta;
    }

    protected void setAtoms() {
        atoms = new AtomList(atom1.molecularSystem);
        atoms.add(atom1);
        atoms.add(atom2);
    }


    public double evaluate() {
        double dis;
        double d, d2, twoD, invDplus, invDplus2;
        double deDd = 0;
        double deDx;
        double deDy;
        double deDz;
        double energy = 0;

        if (frozen()) return 0;
        distance.update();
        dis = distance.distance();
        d = dis - target;
        energy = weight * d * d;
        deDd = 2 * weight * d;
        deDx = deDd * distance.dDistanceDx();
        deDy = deDd * distance.dDistanceDy();
        deDz = deDd * distance.dDistanceDz();
        if (!atom1.frozen()) {
            atom1.addToFx(-1 * deDx); // force = -derivative
            atom1.addToFy(-1 * deDy);
            atom1.addToFz(-1 * deDz);
        }
        if (!atom2.frozen()) {
            atom2.addToFx(deDx);
            atom2.addToFy(deDy);
            atom2.addToFz(deDz);
        }
        return energy;
    }

    public void scaleWeight(double factor) {
        weight *= factor;
    }
}