/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.torsionConstraints;

import meshi.PDB.PdbLine;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneParameters;
import meshi.geometry.Angle;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.QuickAndDirtyTorsion;
import meshi.geometry.Torsion;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.util.UpdateableException;
import meshi.util.Utils;

import java.util.Map;

public class TorsionConstraintElement extends EnergyElement {
    /**
     * The variable BUFF allow a buffer zone around the target, where the energy is flat. Outside
     * this zone the energy rise quadratically. BUFF also control the size of the zone where the
     * two parabolas are fused into derivability at target+PI.
     */
    protected final double BUFF = 0.12;
    protected Atom atom1, atom2, atom3, atom4;
    protected double force, force2, forceSmooth, constSmooth;
    protected double target, std;
    protected Torsion torsion;
    protected double weight;

    protected Angle angle1,angle2;
    protected double[][] angle1Derivatives,angle2Derivatives;
    protected int[][] atomAngleIndices = new int[5][3]; //0 unused; [atom1234][angle12]
    protected double anglesDerivatives[][] = new double[5][3];//0 unused; [atom1234][xyz]


    public TorsionConstraintElement(Torsion torsion, Parameters parameters, double weight) {
        this.weight = weight;
        atom1 = torsion.atom1;
        atom2 = torsion.atom2;
        atom3 = torsion.atom3;
        atom4 = torsion.atom4;
        setAtoms();
        this.torsion = torsion;
        target = ((TorsionConstraintsParameters) parameters).target;
        std  = ((TorsionConstraintsParameters) parameters).std;
    }

    public void setAtoms() {
        atoms = new AtomList(atom1.molecularSystem);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
        atoms.add(atom4);
    }

    private double debugDx1; private double getDebugDx1() {return debugDx1;}
    private double debugDy1; private double getDebugDy1() {return debugDy1;}
    private double debugDx2; private double getDebugDx2() {return debugDx2;}
    public double evaluate() {
        double sinTor = Math.sin(torsion.torsion());
        double cosTor = Math.cos(torsion.torsion());
        double sinTarget = Math.sin(target);
        double cosTarget = Math.cos(target);
        double sinD = sinTor - sinTarget;
        double cosD = cosTor - cosTarget;
//        double gSin = Math.exp(-0.5*(1 /(1 + std)) * sinD * sinD);
//        double gCos = Math.exp(-0.5*(1 /(1 + std)) * cosD * cosD);
//        double torsionEnergy = weight * -1 * gCos*gSin;
//        double torsionDe = -weight * (gCos * gSin * -2 * 0.5 * 1 /(1 + std) * sinD * -1/Math.tan(torsion.torsion()) + gSin * gCos * -2 * 0.5 * 1 /(1 + std) * cosD);
//------------------------
        double torsionEnergy = weight * (sinD*sinD+cosD*cosD) - weight;
        double minusInvTan = -1/Math.tan(torsion.torsion());
        double torsionDe = weight * (2*sinD * minusInvTan + 2*cosD);
//        if ((torsion.getTorsionResNum() == 43) & (Math.abs(Math.tan(torsion.torsion()))< 0.1)) {
//            System.out.println(torsion.name() + "-" + torsion.getTorsionResNum() + " " + 180 * torsion.torsion() / Math.PI + " " +
//                    180 * target / Math.PI + " " + std + " " + (torsion.torsion() - target) + " " + Math.sqrt(sinD * sinD + cosD * cosD) + " " + torsionDe);
//        }
//        if (Math.abs(Math.tan(torsion.torsion()))< 0.00001){
//            System.out.println(torsion.name()+"-"+torsion.getTorsionResNum()+" "+180*torsion.torsion()/Math.PI+" "+
//                    180*target/Math.PI+" "+std+" "+(torsion.torsion()-target)+" "+Math.sqrt(sinD*sinD+cosD*cosD)+" "+torsionDe);
//            throw new RuntimeException("This is a problem");
//        }


        Angle angle1 = torsion.angle1();
        double sinA1 = angle1.sinAngle();
        if (Math.abs(sinA1) < 0.01)
            Utils.println("Unstable angle: "+angle1);
        double minusInvTan1 = -1/Math.tan(angle1.angle());
        double angle1Energy = 1 - Math.exp(-10*sinA1*sinA1*sinA1*sinA1);
        double angle1De    =  40 * (1 - angle1Energy) * sinA1 * sinA1 * sinA1 * minusInvTan1;
//        angle1Energy = 1;
//        angle1De = 0;

        Angle angle2 = torsion.angle2();
        double sinA2 = angle2.sinAngle();
        if (Math.abs(sinA2) < 0.01)
            Utils.println("Unstable angle: "+angle2);
        double minusInvTan2 = -1/Math.tan(angle2.angle());
        double angle2Energy = 1 - Math.exp(-10*sinA2*sinA2*sinA2*sinA2*sinA2*sinA2);
        double ang1e2De    =  40 * (1 - angle2Energy) * sinA2 * sinA2 * sinA2 * sinA2 * sinA2 * minusInvTan2;
//        angle2Energy = 1;
//        ang1e2De = 0;

        double energy = torsionEnergy * angle1Energy * angle2Energy;

        double t_a1_dA2 = torsionEnergy * angle1Energy * ang1e2De;
        double t_dA1_a2 = torsionEnergy * angle1De * angle2Energy;
        double dT_a1_a2 = torsionDe * angle1Energy * angle2Energy;

        double atom1Dx = t_dA1_a2 * angle1.dCosAngleDx1() + dT_a1_a2 * torsion.dCosTorsionDx1(); //Atom1 has no affect on Angle 2
        double atom1Dy = t_dA1_a2 * angle1.dCosAngleDy1() + dT_a1_a2 * torsion.dCosTorsionDy1();
        double atom1Dz = t_dA1_a2 * angle1.dCosAngleDz1() + dT_a1_a2 * torsion.dCosTorsionDz1();

        double atom2Dx = t_a1_dA2 * angle2.dCosAngleDx1() + t_dA1_a2 * angle1.dCosAngleDx2() + dT_a1_a2 * torsion.dCosTorsionDx2();
        debugDx2 = atom2Dx;
        double atom2Dy = t_a1_dA2 * angle2.dCosAngleDy1() + t_dA1_a2 * angle1.dCosAngleDy2() + dT_a1_a2 * torsion.dCosTorsionDy2();
        double atom2Dz = t_a1_dA2 * angle2.dCosAngleDz1() + t_dA1_a2 * angle1.dCosAngleDz2() + dT_a1_a2 * torsion.dCosTorsionDz2();

        double atom3Dx = t_a1_dA2 * angle2.dCosAngleDx2() + t_dA1_a2 * angle1.dCosAngleDx3() + dT_a1_a2 * torsion.dCosTorsionDx3();
        double atom3Dy = t_a1_dA2 * angle2.dCosAngleDy2() + t_dA1_a2 * angle1.dCosAngleDy3() + dT_a1_a2 * torsion.dCosTorsionDy3();
        double atom3Dz = t_a1_dA2 * angle2.dCosAngleDz2() + t_dA1_a2 * angle1.dCosAngleDz3() + dT_a1_a2 * torsion.dCosTorsionDz3();

        double atom4Dx = t_a1_dA2 * angle2.dCosAngleDx3() + dT_a1_a2 * torsion.dCosTorsionDx4();
        double atom4Dy = t_a1_dA2 * angle2.dCosAngleDy3() + dT_a1_a2 * torsion.dCosTorsionDy4();
        double atom4Dz = t_a1_dA2 * angle2.dCosAngleDz3() + dT_a1_a2 * torsion.dCosTorsionDz4();

        atom1.addToFx(-atom1Dx);
        atom1.addToFy(-atom1Dy);
        atom1.addToFz(-atom1Dz);

        atom2.addToFx(-atom2Dx);
        atom2.addToFy(-atom2Dy);
        atom2.addToFz(-atom2Dz);

        atom3.addToFx(-atom3Dx);
        atom3.addToFy(-atom3Dy);
        atom3.addToFz(-atom3Dz);

        atom4.addToFx(-atom4Dx);
        atom4.addToFy(-atom4Dy);
        atom4.addToFz(-atom4Dz);

        return energy;
    }


    public void scaleWeight(double factor) {
        weight *= factor;
    }
//-------------------------------------------------------------------------------------------------------------------------------------
    public static void main(String[] args) throws UpdateableException {
        double epsilon = 0.0000000001;





        PdbLine pdbLine1 = new PdbLine("ATOM   2156  CA  SER   242      61.081   0.362   0.280  1.00  0.00");
        PdbLine pdbLine2 = new PdbLine("ATOM   2160  C   SER   242      59.972  -0.632  -0.048  1.00  0.00");
        PdbLine pdbLine3 = new PdbLine("ATOM   2164  N   THR   243      58.838  -0.507   0.634  1.00  0.00");
        PdbLine pdbLine4 = new PdbLine("ATOM   2163  CA  THR   243      57.641  -1.247   0.253  1.00  0.00");

        //Torsion torsion = new Torsion()

        MolecularSystem ms = new MolecularSystem();
        Atom atom1 = new Atom(pdbLine1, ms);
        Atom atom2 = new Atom(pdbLine2, ms);
        Atom atom3 = new Atom(pdbLine3, ms);
        Atom atom4 = new Atom(pdbLine4, ms);
        ms.terminator().reset();
        DistanceMatrix distanceMatrix = new DistanceMatrix(ms);
        Angle angle1 = new Angle(atom1, atom2, atom3, distanceMatrix);
        Angle angle2 = new Angle(atom2, atom3, atom4, distanceMatrix);
        Torsion torsion = new QuickAndDirtyTorsion(angle1, angle2, distanceMatrix);

        OutOfPlaneParameters parameters= new OutOfPlaneParameters(null, null, null, null, 80, 1);
        updateTest(angle1, angle2, torsion, distanceMatrix);
        TorsionConstraintElement element = new TorsionConstraintElement(torsion, parameters, 1);

        updateTest(angle1, angle2, torsion, distanceMatrix);

        element.evaluate();
        double c1 = element.evaluate();
        double a = element.getDebugDx2();
        atom2.addToX(epsilon);
        updateTest(angle1, angle2, torsion, distanceMatrix);
        element.evaluate();
        double c2 = element.evaluate();
        System.out.println(torsion.torsion()+" "+angle1.angle()+" "+angle2.angle()+" "+c1+" "+a+" "+(c2-c1)/epsilon);
    }

    private static void updateTest(Angle angle1, Angle angle2, Torsion torsion, DistanceMatrix distanceMatrix) throws UpdateableException{
        distanceMatrix.update();
        angle1.update();
        angle2.update();
        torsion.update();
    }

    public void report() {
        System.out.println(torsion.name()+"-"+torsion.getTorsionResNum()+" "+180*torsion.torsion()/Math.PI+" "+180*target/Math.PI+" "+std+" "+(torsion.torsion()-target)+" "+evaluate()/weight);
    }

}
