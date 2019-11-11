/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.plane;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.geometry.Angle;
import meshi.geometry.Torsion;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;

//---------------------------------------------------------------------
public class PlaneEnergyElement extends EnergyElement {
    public final double BETA = 0.05;
    protected Atom atom1, atom2, atom3, atom4;
    protected double force, force2;
    protected PlaneParameters.Isomer isomer;
    protected Torsion torsion;
    protected double weight;
    protected Angle angle1,angle2;
    protected double[][] angle1Derivatives,angle2Derivatives;
    protected int[][] atomAngleIndices = new int[5][3]; //0 unused; [atom1234][angle12]
    protected double anglesDerivatives[][] = new double[5][3];//0 unused; [atom1234][xyz]
    private boolean on;

    public PlaneEnergyElement(Torsion torsion, Parameters parameters, double weight){
        this.weight = weight;
        atom1 = torsion.atom1;
        atom2 = torsion.atom2;
        atom3 = torsion.atom3;
        atom4 = torsion.atom4;
        this.torsion = torsion;
        setAtoms();
        updateFrozen();
        on = true;
        isomer = ((PlaneParameters) parameters).isomer;

        force = ((PlaneParameters) parameters).force * weight;
        force2 = ((PlaneParameters) parameters).force2 * weight;
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

    public void off() {on = false;}
    public void setAtoms() {
        atoms = new AtomList(atom1.molecularSystem);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
        atoms.add(atom4);
    }

    public PlaneParameters.Isomer getIsomer() {
        return isomer;
    }
    public Torsion getTorsion() {
        return torsion;
    }



    public void scaleWeight(double factor) {
        force *= factor;
        force2 *= factor;
    }

    public double evaluate() {
        double dEdCosTorsion;
        double cosTorsion;
        double ctpo; // Cos Torsion Plus One
        double ctpo2; // (Cos Torsion Plus One)^2
        double invCtpo2pb;
        double omct; // One Minus Cos Torsion
        double omct2;
        double invOmct2pb;
        double energy, energy1 = 0;
        double anglesEnergy;


        if (frozen()) return 0;
        if (!on) return 0;

        cosTorsion = torsion.cosTorsion();
         if (isomer == PlaneParameters.Isomer.TRANS) {
            ctpo = cosTorsion + 1; // Cos Torsion Plus One
            ctpo2 = ctpo * ctpo;
            invCtpo2pb = 1 / (ctpo2 + BETA);
            energy1 = force * (ctpo2 * invCtpo2pb);
            dEdCosTorsion = force2 * ctpo * (invCtpo2pb - ctpo2 * invCtpo2pb * invCtpo2pb);

         } else if (isomer == PlaneParameters.Isomer.CIS) {
            omct = 1 - cosTorsion; // One Minus Cos Torsion
            omct2 = omct * omct;
            invOmct2pb = 1 / (omct2 + BETA);
            energy1 = force * (omct2 + omct2 * invOmct2pb);
            dEdCosTorsion = -1 * force2 * omct * (1 + invOmct2pb - omct2 * invOmct2pb * invOmct2pb);
         } else {
            energy1 = force * (1 - cosTorsion * cosTorsion);
            dEdCosTorsion = -1 * force2 * cosTorsion;
         }

        if ((!(energy1> 0)) && (!(energy1 == 0)) &&(!(energy1<0)) )
            throw new RuntimeException("weird angle energy element 2 "+energy1+" "+cosTorsion+" "+isomer+"\n"+torsion+"\n"+angle1+"\n"+angle2);

        anglesEnergy = anglesEnergy(angle1,angle2, atomAngleIndices,angle1Derivatives,angle2Derivatives, anglesDerivatives);
        energy = energy1*anglesEnergy;
    //    if (energy1 > 100000)
    //        throw new RuntimeException(torsion+"\n"+energy1+" "+anglesEnergy+" "+energy+"  ******************\n");
        if (!atom1.frozen()) {
            atom1.addToFx(-1 * dEdCosTorsion * torsion.dCosTorsionDx1()* anglesEnergy);
            atom1.addToFy(-1 * dEdCosTorsion * torsion.dCosTorsionDy1()* anglesEnergy);
            atom1.addToFz(-1 * dEdCosTorsion * torsion.dCosTorsionDz1()* anglesEnergy);
            atom1.addToFx(-1 * energy1 * anglesDerivatives[1][0]);
            atom1.addToFy(-1 * energy1 * anglesDerivatives[1][1]);
            atom1.addToFz(-1 * energy1 * anglesDerivatives[1][2]);
        }
        if (!atom2.frozen()) {
            atom2.addToFx(-1 * dEdCosTorsion * torsion.dCosTorsionDx2()* anglesEnergy);
            atom2.addToFy(-1 * dEdCosTorsion * torsion.dCosTorsionDy2()* anglesEnergy);
            atom2.addToFz(-1 * dEdCosTorsion * torsion.dCosTorsionDz2()* anglesEnergy);
            atom2.addToFx(-1 * energy1 * anglesDerivatives[2][0]);
            atom2.addToFy(-1 * energy1 * anglesDerivatives[2][1]);
            atom2.addToFz(-1 * energy1 * anglesDerivatives[2][2]);
        }
        if (!atom3.frozen()) {
            atom3.addToFx(-1 * dEdCosTorsion * torsion.dCosTorsionDx3()* anglesEnergy);
            atom3.addToFy(-1 * dEdCosTorsion * torsion.dCosTorsionDy3()* anglesEnergy);
            atom3.addToFz(-1 * dEdCosTorsion * torsion.dCosTorsionDz3()* anglesEnergy);
            atom3.addToFx(-1 * energy1 * anglesDerivatives[3][0]);
            atom3.addToFy(-1 * energy1 * anglesDerivatives[3][1]);
            atom3.addToFz(-1 * energy1 * anglesDerivatives[3][2]);
        }
        if (!atom4.frozen()) {
            atom4.addToFx(-1 * dEdCosTorsion * torsion.dCosTorsionDx4()* anglesEnergy);
            atom4.addToFy(-1 * dEdCosTorsion * torsion.dCosTorsionDy4()* anglesEnergy);
            atom4.addToFz(-1 * dEdCosTorsion * torsion.dCosTorsionDz4()* anglesEnergy);
            atom4.addToFx(-1 * energy1 * anglesDerivatives[4][0]);
            atom4.addToFy(-1 * energy1 * anglesDerivatives[4][1]);
            atom4.addToFz(-1 * energy1 * anglesDerivatives[4][2]);
        }
        return energy;
    }

    public String toString() {
        String prop = ((torsion.proper()) ? "proper" : "improper");
        return ("PlaneEnergyElement " + torsion.name() + " " + isomer + " " + prop + "\n" +
                "force = " + dFormatSrt.f(force) + " torsion = " +
                dFormatSrt.f(Angle.rad2deg(torsion.torsion())) + " energy = " + dFormatSrt.f(evaluate()) + "\n" +
                atom1.verbose(1) + "\n" + atom2.verbose(1) + "\n" + atom3.verbose(1) + "\n" + atom4.verbose(1)) + "\n" + torsion;
    }

    // for numerical stability we want the derivative of the energy to be zero at PI
    public static void angleEnergy(Angle angle,double[] out) {
        double epsilon = 0.1;
        double alpha = 100;
        double angleMinusPI = angle.angle()-Math.PI;
        double angleMinusPiSquare = alpha*angleMinusPI*angleMinusPI;
        double angleMinusPi2 = 2*alpha*angleMinusPI;
        double angleMinusPiSquarePlusEpsilon = angleMinusPiSquare+epsilon;
        double e = angleMinusPiSquare/angleMinusPiSquarePlusEpsilon;
        double de = (angleMinusPi2/angleMinusPiSquarePlusEpsilon)*(1-e);
        out[0] = e;
        out[1] = de;
    }

    public static double anglesEnergy( Angle angle1, Angle angle2, int[][] atomAngleIndices,
                                     double[][] angle1Derivatives,
                                     double[][] angle2Derivatives,
                                     double[][] anglesDerivatives) {
        double e1,e2,d1,d2;
        double[] angle1energy = new double[2];  //energy; derivative;
        double[] angle2energy = new double[2];
        double anglesEnergy;

        int atom1Angle1 = atomAngleIndices[1][1];
        int atom1Angle2 = atomAngleIndices[1][2];
        int atom2Angle1 = atomAngleIndices[2][1];
        int atom2Angle2 = atomAngleIndices[2][2];
        int atom3Angle1 = atomAngleIndices[3][1];
        int atom3Angle2 = atomAngleIndices[3][2];
        int atom4Angle1 = atomAngleIndices[4][1];
        int atom4Angle2 = atomAngleIndices[4][2];
        angleEnergy(angle1,angle1energy);
        angleEnergy(angle2,angle2energy);
        e1 = angle1energy[0];
        e2 = angle2energy[0];
        d1 = angle1energy[1];
        d2 = angle2energy[1];
        anglesEnergy = e1*e2;
        if ((!(anglesEnergy> 0)) && (!(anglesEnergy == 0)) &&(!(anglesEnergy<0)) )
            throw new RuntimeException("weird angle energy element 1"+anglesEnergy+"\n"+angle1+"\n"+angle2);

        anglesDerivatives[1][0] = e1*d2*angle2Derivatives[atom1Angle2][0];
        anglesDerivatives[1][1] = e1*d2*angle2Derivatives[atom1Angle2][1];
        anglesDerivatives[1][2] = e1*d2*angle2Derivatives[atom1Angle2][2];
        anglesDerivatives[1][0] += e2*d1*angle1Derivatives[atom1Angle1][0];
        anglesDerivatives[1][1] += e2*d1*angle1Derivatives[atom1Angle1][1];
        anglesDerivatives[1][2] += e2*d1*angle1Derivatives[atom1Angle1][2];

        anglesDerivatives[2][0] = e1*d2*angle2Derivatives[atom2Angle2][0];
        anglesDerivatives[2][1] = e1*d2*angle2Derivatives[atom2Angle2][1];
        anglesDerivatives[2][2] = e1*d2*angle2Derivatives[atom2Angle2][2];
        anglesDerivatives[2][0] += e2*d1*angle1Derivatives[atom2Angle1][0];
        anglesDerivatives[2][1] += e2*d1*angle1Derivatives[atom2Angle1][1];
        anglesDerivatives[2][2] += e2*d1*angle1Derivatives[atom2Angle1][2];

        anglesDerivatives[3][0] = e1*d2*angle2Derivatives[atom3Angle2][0];
        anglesDerivatives[3][1] = e1*d2*angle2Derivatives[atom3Angle2][1];
        anglesDerivatives[3][2] = e1*d2*angle2Derivatives[atom3Angle2][2];
        anglesDerivatives[3][0] += e2*d1*angle1Derivatives[atom3Angle1][0];
        anglesDerivatives[3][1] += e2*d1*angle1Derivatives[atom3Angle1][1];
        anglesDerivatives[3][2] += e2*d1*angle1Derivatives[atom3Angle1][2];

        anglesDerivatives[4][0] = e1*d2*angle2Derivatives[atom4Angle2][0];
        anglesDerivatives[4][1] = e1*d2*angle2Derivatives[atom4Angle2][1];
        anglesDerivatives[4][2] = e1*d2*angle2Derivatives[atom4Angle2][2];
        anglesDerivatives[4][0] += e2*d1*angle1Derivatives[atom4Angle1][0];
        anglesDerivatives[4][1] += e2*d1*angle1Derivatives[atom4Angle1][1];
        anglesDerivatives[4][2] += e2*d1*angle1Derivatives[atom4Angle1][2];
        return anglesEnergy;
    }

    public Atom getAtom1() {
        return atom1;
    }

    public Atom getAtom2() {
        return atom2;
    }

    public Atom getAtom3() {
        return atom3;
    }

    public Atom getAtom4() {
        return atom4;
    }

}
