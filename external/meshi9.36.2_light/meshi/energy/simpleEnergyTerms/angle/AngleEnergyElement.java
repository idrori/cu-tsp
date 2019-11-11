/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

//-----------------------------------------------------------------------------------------------
package meshi.energy.simpleEnergyTerms.angle;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.geometry.Angle;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;

public class AngleEnergyElement extends EnergyElement  {
    protected Atom atom1, atom2, atom3;
    protected int number1, number2, number3;
    protected Angle angle;
    public double target, force, force2;
    double weight;

    public AngleEnergyElement(Angle angle, Parameters parameters, double weight) {
        this.weight = weight;
        atom1 = angle.atom1;
        atom2 = angle.atom2;
        atom3 = angle.atom3;
        this.angle = angle;
        target = ((AngleParameters) parameters).target;
        if ((target < Angle.TOO_LOW) || (target > Angle.TOO_HIGH))
            throw new RuntimeException("The target angle " + target + " is too close to 0 or PI." +
                    "currently MESHI can not handle those values.");
        force = ((AngleParameters) parameters).force * weight;
        force2 = ((AngleParameters) parameters).force2 * weight;
        setAtoms();
        updateFrozen();
    }

    protected void setAtoms() {
        atoms = new AtomList(atom1.molecularSystem);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
    }

    /**
     * Angle energy calculation and atom forces updating.
     */
    public double evaluate(){// throws EvaluationException{
        double energy = 0;
        double angleValue;

        angleValue = angle.angle();
//        if (angle.dangerous())
//            throw new InstableAngleException("\n Angle Energy element of\n "+angle+"\n"+
//                    180*angleValue/Math.PI+"degrees  is either too small or too large for stable numerical processing.\n");
            if (frozen()) return 0;
            energy = calcEnergy(atom1,atom2,atom3,angleValue,target,force);
//            if (angleValue <= Angle.TOO_LOW)
//                energy += calcEnergy(atom1,atom2,atom3,angleValue, Angle.TOO_LOW,100*force);
//            else if (angleValue >= Angle.TOO_HIGH) {
//                System.out.println("\n*********\n"+angle+"\n");
//                energy += calcEnergy(atom1,atom2,atom3,angleValue, Angle.TOO_HIGH,100*force);
//            }
        return energy;
    }
    private double calcEnergy(Atom atom1, Atom atom2, Atom atom3,
                              double angleValue, double target, double force ){
        double d;
        double deDd;
        double energy;
            d = angleValue - target;
            energy = d * d * force;
            deDd = -1 * d * force2; // force = -derivative
            if (!atom1.frozen()) {
                atom1.addToFx(deDd * angle.dangleDx1());
                atom1.addToFy(deDd * angle.dangleDy1());
                atom1.addToFz(deDd * angle.dangleDz1());
            }
            if (!atom2.frozen()) {
                atom2.addToFx(deDd * angle.dangleDx2());
                atom2.addToFy(deDd * angle.dangleDy2());
                atom2.addToFz(deDd * angle.dangleDz2());
            }
            if (!atom3.frozen()) {
                atom3.addToFx(deDd * angle.dangleDx3());
                atom3.addToFy(deDd * angle.dangleDy3());
                atom3.addToFz(deDd * angle.dangleDz3());
            }

        return energy;
    }


    public String toString() {
        double energy;
  //      try {
            energy = evaluate();
//        }
//        catch (EvaluationException ex) {
//            System.out.println(ex);
//            energy = 9999999;
//        }
        return ("AngleEnergyElement target = " + dFormatSrt.f(Angle.rad2deg(target)) + " force = " + dFormatSrt.f(force) + " angle = " +
                dFormatSrt.f(Angle.rad2deg(angle.angle())) + " energy = " + dFormatSrt.f(energy) + "\n" +
                atom1.verbose(1) + "\n" + atom2.verbose(1) + "\n" + atom3.verbose(1));
    }
}
