/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.outOfPlane;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.plane.PlaneEnergyElement;
import meshi.geometry.Angle;
import meshi.geometry.Torsion;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;

public class OutOfPlaneEnergyElement extends EnergyElement {
    /**
     * The variable BUFF allow a buffer zone around the target, where the energy is flat. Outside
     * this zone the energy rise quadratically. BUFF also control the size of the zone where the
     * two parabolas are fused into derivability at target+PI.
     */
    public static final double BUFF = 0.12;
    protected Atom atom1, atom2, atom3, atom4;
    protected double target, force, force2, forceSmooth, constSmooth;
    protected double targetMinusBUFF, targetPlusBUFF, targetPlusPI, targetPlusPIMinusBUFF, targetPlusPIPlusBUFF, targetPlus2PIMinusBUFF;
    protected Torsion torsion;
    protected double weight;

    protected Angle angle1,angle2;
    protected double[][] angle1Derivatives,angle2Derivatives;
    protected int[][] atomAngleIndices = new int[5][3]; //0 unused; [atom1234][angle12]
    protected double anglesDerivatives[][] = new double[5][3];//0 unused; [atom1234][xyz]


    public OutOfPlaneEnergyElement(Torsion torsion, Parameters parameters, double weight) {
        this.weight = weight;
        atom1 = torsion.atom1;
        atom2 = torsion.atom2;
        atom3 = torsion.atom3;
        atom4 = torsion.atom4;
        setAtoms();
        this.torsion = torsion;
        target = ((OutOfPlaneParameters) parameters).target;
        if (((target - BUFF) < -Math.PI) || ((target + BUFF) > 0))
            throw new RuntimeException("Currently, the functional format of OutOfPlaneEnergyElement can only\n" +
                    "support targets that are smaller than -BUFF radians and larger than -PI+buff. If this pose a \n" +
                    "problem adjust the functional form in the evaluate() method\n"+
            "target = "+target+"\n"+
            "BUFF = "+BUFF);
        force = ((OutOfPlaneParameters) parameters).force * weight;
        force2 = ((OutOfPlaneParameters) parameters).force2 * weight;
        targetMinusBUFF = target - BUFF;
        targetPlusBUFF = target + BUFF;
        targetPlusPI = target + Math.PI;
        targetPlusPIMinusBUFF = target + Math.PI - BUFF;
        targetPlusPIPlusBUFF = target + Math.PI + BUFF;
        targetPlus2PIMinusBUFF = target + 2 * Math.PI - BUFF;
        forceSmooth = -force * (Math.PI - 2 * BUFF) / BUFF;
        constSmooth = force * (Math.PI - 2 * BUFF) * (Math.PI - BUFF);
        updateFrozen();

        angle1 = torsion.angle1();
        angle2 = torsion.angle2();
        angle1Derivatives = torsion.angle1Derivatives();
        angle2Derivatives = torsion.angle2Derivatives();

        atomAngleIndices[1][1] = torsion.getAtom1Angle1();
        atomAngleIndices[1][2] = torsion.getAtom1Angle2();
        atomAngleIndices[2][1] = torsion.getAtom2Angle1();
        atomAngleIndices[2][2] = torsion.getAtom2Angle2();
        atomAngleIndices[3][1] = torsion.getAtom3Angle1();
        atomAngleIndices[3][2] = torsion.getAtom3Angle2();
        atomAngleIndices[4][1] = torsion.getAtom4Angle1();
        atomAngleIndices[4][2] = torsion.getAtom4Angle2();

    }

    public void setAtoms() {
        atoms = new AtomList(atom1.molecularSystem);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
        atoms.add(atom4);
    }

    public double evaluate() {
        double dEdTorsion, e;
        double torsionValue, dTorsion;
        double energy1 = 0;
        double energy;
        double anglesEnergy;

        if (frozen()) return 0;
        torsionValue = torsion.torsion();
        if ((torsionValue > targetMinusBUFF) && (torsionValue < targetPlusBUFF)) {
            dEdTorsion = 0.0;
        } else if (torsionValue <= targetMinusBUFF) {
            dTorsion = torsionValue - targetMinusBUFF;
            energy1 += force * dTorsion * dTorsion;
            dEdTorsion = force2 * dTorsion;
        } else if ((torsionValue >= targetPlusBUFF) && (torsionValue < targetPlusPIMinusBUFF)) {
            dTorsion = torsionValue - targetPlusBUFF;
            energy1 += force * dTorsion * dTorsion;
            dEdTorsion = force2 * dTorsion;
        } else if (torsionValue > targetPlusPIPlusBUFF) {
            dTorsion = torsionValue - targetPlus2PIMinusBUFF;
            energy1 += force * dTorsion * dTorsion;
            dEdTorsion = force2 * dTorsion;
        } else if ((torsionValue >= targetPlusPIMinusBUFF) && (torsionValue <= targetPlusPIPlusBUFF)) {
            dTorsion = torsionValue - targetPlusPI;
            energy1 += forceSmooth * dTorsion * dTorsion + constSmooth;
            dEdTorsion = 2 * forceSmooth * dTorsion;
        } else {
            throw new RuntimeException("If the run got here, there is a bug in the code.\n" +
                    "torsion is" + torsion + "\n" +
                    "torsion Value = " + torsionValue + "\n" +
                    "targetMinusBUFF & targetPlusBUFF " + targetMinusBUFF + "  " + targetPlusBUFF + "\n" +
                    "targetPlusPIMinusBUFF " + targetPlusPIMinusBUFF);
        }

        anglesEnergy = PlaneEnergyElement.anglesEnergy(angle1, angle2, atomAngleIndices,
                                                       angle1Derivatives, angle2Derivatives,
                                                       anglesDerivatives);
        energy = energy1*anglesEnergy;
        if (!atom1.frozen()) {
            atom1.addToFx(-1 * dEdTorsion * torsion.dTorsionDx1()* anglesEnergy);
            atom1.addToFy(-1 * dEdTorsion * torsion.dTorsionDy1()* anglesEnergy);
            atom1.addToFz(-1 * dEdTorsion * torsion.dTorsionDz1()* anglesEnergy);
            atom1.addToFx(-1 * energy1 * anglesDerivatives[1][0]);
            atom1.addToFy(-1 * energy1 * anglesDerivatives[1][1]);
            atom1.addToFz(-1 * energy1 * anglesDerivatives[1][2]);
        }
        if (!atom2.frozen()) {
            atom2.addToFx(-1 * dEdTorsion * torsion.dTorsionDx2()* anglesEnergy);
            atom2.addToFy(-1 * dEdTorsion * torsion.dTorsionDy2()* anglesEnergy);
            atom2.addToFz(-1 * dEdTorsion * torsion.dTorsionDz2()* anglesEnergy);
            atom2.addToFx(-1 * energy1 * anglesDerivatives[2][0]);
            atom2.addToFy(-1 * energy1 * anglesDerivatives[2][1]);
            atom2.addToFz(-1 * energy1 * anglesDerivatives[2][2]);
        }
        if (!atom3.frozen()) {
            atom3.addToFx(-1 * dEdTorsion * torsion.dTorsionDx3()* anglesEnergy);
            atom3.addToFy(-1 * dEdTorsion * torsion.dTorsionDy3()* anglesEnergy);
            atom3.addToFz(-1 * dEdTorsion * torsion.dTorsionDz3()* anglesEnergy);
            atom3.addToFx(-1 * energy1 * anglesDerivatives[3][0]);
            atom3.addToFy(-1 * energy1 * anglesDerivatives[3][1]);
            atom3.addToFz(-1 * energy1 * anglesDerivatives[3][2]);
        }
        if (!atom4.frozen()) {
            atom4.addToFx(-1 * dEdTorsion * torsion.dTorsionDx4()* anglesEnergy);
            atom4.addToFy(-1 * dEdTorsion * torsion.dTorsionDy4()* anglesEnergy);
            atom4.addToFz(-1 * dEdTorsion * torsion.dTorsionDz4()* anglesEnergy);
            atom4.addToFx(-1 * energy1 * anglesDerivatives[4][0]);
            atom4.addToFy(-1 * energy1 * anglesDerivatives[4][1]);
            atom4.addToFz(-1 * energy1 * anglesDerivatives[4][2]);
        }
        return energy;
    }

    public String toString() {
        return ("OutOfPlaneEnergyElement " + torsion.name() + " target = " + Angle.rad2deg(target) +
                " force = " + dFormatSrt.f(force) + " torsion = " +
                dFormatSrt.f(Angle.rad2deg(torsion.torsion())) + " energy = " + dFormatSrt.f(evaluate()) + "\n" +
                atom1.verbose(1) + "\n" + atom2.verbose(1) + "\n" + atom3.verbose(1) + "\n" + atom4.verbose(1));
    }
}
