/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.TotalEnergy;
import meshi.energy.EnergyElement;
import meshi.geometry.DistanceLists;
import meshi.parameters.AtomType;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.DistanceList;
import meshi.geometry.Distance;
import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.MolecularSystem;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 29/06/2009
 * Time: 10:52:58
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZSolvate extends AbstractEnergy {

    private SolvateEnergy solvateTerm;
    public static final double[] mean = new double[AtomType.values().length];
    public static final double[] std = new double[AtomType.values().length];
    public static final double EPSILON = 0.1;


    private static DistanceMatrix distanceMatrix;

    public double factor, e3;
    private MolecularSystem molecularSystem;

    public double weight() {
        return weight;
    }

    public CooperativeZSolvate(DistanceMatrix distanceMatrix, CooperativeZSolvateParameters parameters,
                               SolvateEnergy solvateTerm, EnergyInfoElement info) {
        super(toArray(distanceMatrix), info);
        if (weight != 0) {
            this.solvateTerm = solvateTerm;
            if (solvateTerm == null) throw new RuntimeException("Cooperative Solvate term could not be created " +
                    "before the non-cooperative term is created.");
            comment = "CooperativeZSolvate";
            for (AtomType type : AtomType.values()) {
                mean[type.ordinal()] = parameters.mean[type.ordinal()];
                std[type.ordinal()] = parameters.std[type.ordinal()];
            }
            this.distanceMatrix = distanceMatrix;
            molecularSystem = distanceMatrix.molecularSystem;
        }
    }

    public EnergyInfoElement evaluate() {

        return evaluate(false);
    }

    public void evaluateAtoms() {
        evaluate(true);
    }

    public EnergyInfoElement evaluate(boolean evaluateAtoms) {
        DistanceLists nonBondedList = distanceMatrix.nonBondedList();
        double energy, e, e3, e2, eMinusEpsilon, eMinusEpsilon2;
        double derivative, derivative1, derivative2;
        AtomType type1, type2, type;
        e = 0;
        for (AtomCore atom : molecularSystem) {
            type = atom.type();
            e += zScore(solvateTerm.atomEnergies[atom.number], mean[type.ordinal()], std[type.ordinal()]);
        }
        if (e > 0) {
            info.setValue(0);
            return info;
        }
        e2 = e * e;
        e3 = e2 * e;
        eMinusEpsilon = e - EPSILON;
        eMinusEpsilon2 = eMinusEpsilon * eMinusEpsilon;
        energy = weight * e3 / eMinusEpsilon;
        derivative = weight * (3 * e2 / eMinusEpsilon - e3 / eMinusEpsilon2);
        if (evaluateAtoms) {
            for (AtomCore atom : molecularSystem) {
                type = atom.type();
                atom.atom.addEnergy(energy * zScore(solvateTerm.atomEnergies[atom.number], mean[type.ordinal()], std[type.ordinal()]) / e);
            }
        }


        double forceX1, forceY1, forceZ1, forceX2, forceY2, forceZ2;
        forceX1 = forceY1 = forceZ1 = forceX2 = forceY2 = forceZ2 = 0;

        int atomNumber1, atomNumber2;
        SolvateDistanceAttribute sigmaValues = null;
        SolvateDistanceAttributeBetweenPolars sigmaValuesBP = null;
        SolvateDistanceAttributeWithNonPolar sigmaValuesWNP = null;
        SolvateDistanceAttributeNonPolarPolar sigmaValuesNPP = null;
        SolvateDistanceAttributePolarNonPolar sigmaValuesPNP = null;
        SolvateDistanceType typeOfSolvateDis;
        AtomCore atom1, atom2;

        for (DistanceList distanceList : nonBondedList) {
            for (Distance dis : distanceList) {
                sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATION_ALL_ATOM_ATTRIBUTE);
                typeOfSolvateDis = sigmaValues.type();
                if (sigmaValues == SolvateEnergy.NOT_PATRICIPATING_DISTANCE) continue;

                atom1 = dis.atom1;
                atom2 = dis.atom2;

                atomNumber1 = atom1.number;
                atomNumber2 = atom2.number;

                type1 = atom1.type();
                type2 = atom2.type();
                if (std[type1.ordinal()] != 0) derivative1 = derivative / std[type1.ordinal()];
                else derivative1 = 0;
                if (std[type2.ordinal()] != 0) derivative2 = derivative / std[type2.ordinal()];
                else derivative2 = 0;
                // if ((std[type1.ordinal()] == 0.0) || (std[type2.ordinal()] == 0.0))
                //       continue;
                //                              throw new RuntimeException("Something is weird in Parameters "+this);
            /*
                        int a = 0;
                        if ((atomNumber1 == 6) || (atomNumber2 == 6))
                        a = 6;

                        if ((atomNumber1 == 4) || (atomNumber2 == 4))
                        a = 4;
            //*/
                switch (typeOfSolvateDis) {
                    case HYDROPHOBIC_POLAR:
                        sigmaValuesNPP = (SolvateDistanceAttributeNonPolarPolar) sigmaValues;
                        forceX1 = -solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dx * derivative2;
                        forceY1 = -solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dy * derivative2;
                        forceZ1 = -solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dz * derivative2;
                        forceX2 = solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dx * derivative2;
                        forceY2 = solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dy * derivative2;
                        forceZ2 = solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dz * derivative2;

                        if (!atom1.status().frozen()) atom1.addForce(forceX1, forceY1, forceZ1);
                        if (!atom2.status().frozen()) atom2.addForce(forceX2, forceY2, forceZ2);
                        break;
                    case POLAR_HYOPHOBIC:
                        sigmaValuesPNP = (SolvateDistanceAttributePolarNonPolar) sigmaValues;
                        forceX1 = -solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dx * derivative1;
                        forceY1 = -solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dy * derivative1;
                        forceZ1 = -solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dz * derivative1;
                        forceX2 = solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dx * derivative1;
                        forceY2 = solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dy * derivative1;
                        forceZ2 = solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dz * derivative1;

                        if (!atom1.status().frozen()) atom1.addForce(forceX1, forceY1, forceZ1);
                        if (!atom2.status().frozen()) atom2.addForce(forceX2, forceY2, forceZ2);
                        break;
                    case TWO_HYDROPHOBIC:
                        sigmaValuesWNP = (SolvateDistanceAttributeWithNonPolar) sigmaValues;
                        // Doing the self derivatives
                        forceX1 = -solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dx * derivative1;
                        forceY1 = -solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dy * derivative1;
                        forceZ1 = -solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dz * derivative1;

                        forceX2 = solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dx * derivative2;
                        forceY2 = solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dy * derivative2;
                        forceZ2 = solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dz * derivative2;

                        // Doing the cross derivatives
                        forceX2 += solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dx * derivative1;
                        forceY2 += solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dy * derivative1;
                        forceZ2 += solvateTerm.dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dz * derivative1;

                        forceX1 -= solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dx * derivative2;
                        forceY1 -= solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dy * derivative2;
                        forceZ1 -= solvateTerm.dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dz * derivative2;

                        if (!atom1.status().frozen()) atom1.addForce(forceX1, forceY1, forceZ1);
                        if (!atom2.status().frozen()) atom2.addForce(forceX2, forceY2, forceZ2);
                        break;
                    case TWO_POLARS:
                        sigmaValuesBP = (SolvateDistanceAttributeBetweenPolars) sigmaValues;
                        if (sigmaValuesBP.sigmHB > 0.0) { // The HB related derivatives
                            double solvFactor = solvateTerm.dSplineDHBC[atomNumber1] * sigmaValuesBP.saltBridgeFactorA1 * derivative1 +
                                    solvateTerm.dSplineDHBC[atomNumber2] * sigmaValuesBP.saltBridgeFactorA2 * derivative2;
                            solvateTerm.solvateHB.findBondByPolars(dis.atom1(), dis.atom2()).applyForcesToAtoms(solvFactor);
                        }
                        break;
                }


                //          if(!atom1.status().frozen() ) atom1.addForce(-factor*forceX1,  -factor*forceY1,  -factor*forceZ1);
                //        if( !atom2.status().frozen() ) atom2.addForce(-factor*forceX2, -factor*forceY2, -factor*forceZ2);


            }
//               */
        /*    Atom atom;
        double diff;
        for (int i = 0 ; i<molecularSystem.size() ; i++) {
            if ( molecularSystem.get(i).atom.nowhere()) continue;
            atom = molecularSystem.get(i).atom;
            if (! atom.frozen()) {
                type = atom.type();
                if (std[type.ordinal()] == 0.0)
                    continue;
                                        //throw new RuntimeException("Something is weird in Parameters "+this);
                //diff     = 2*e/std[type.ordinal()];
                  diff     = 2*e;
                factor    = diff*weight;
                if (factor != 0) {
                    atom.addToFx(-solvateTerm.forceX[i]*factor);
                    atom.addToFy(-solvateTerm.forceY[i]*factor);
                    atom.addToFz(-solvateTerm.forceZ[i]*factor);
                }
            }

        }
        //*/
        }

        info.setValue(energy);
        return info;
    }


    private double zScore(double curE, double mean, double std) {
        // double koeffOutStd = 1;
        if (std == 0)
            return 0;
        //      throw new RuntimeException("Something weird in Parameters "+this);
        /*       if (Math.abs(mean - curE) < std ){
                 koeffOutStd = DEFAULT_KOEFF_INSTD;
             }
             else
                 koeffOutStd = DEFAULT_KOEFF_OUTSTD;
        */

        return (curE - mean) / std;    // TODO fix this

    }


    public void test(TotalEnergy totalEnergy, Atom atom) {
        if (!on) {
            System.out.println("" + this + " is off");
            return;
        }
        System.out.println("Testing " + this + " using " + atom);
        if (atom == null)
            throw new RuntimeException("Cannot test " + this);

        double[][] coordinates = new double[3][];
        coordinates[0] = atom.X();
        coordinates[1] = atom.Y();
        coordinates[2] = atom.Z();
        for (int i = 0; i < 3; i++) {
                totalEnergy.update();

            solvateTerm.evaluate();
            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            double e1 = evaluate().energy();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += EnergyElement.DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
                totalEnergy.update();
            solvateTerm.evaluate();
            double e2 = evaluate().energy();
            double de = e2 - e1;
            double numericalForce = -de / EnergyElement.DX;
            coordinates[i][0] -= EnergyElement.DX;
                totalEnergy.update();

            double diff = Math.abs(analiticalForce - numericalForce);

            if ((2 * diff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + EnergyElement.VERY_SMALL)) > EnergyElement.relativeDiffTolerance) {
                System.out.println("Testing " + this);
                System.out.println("Atom[" + atom.number() + "]." + EnergyElement.XYZ.charAt(i) + " = " + x);
                System.out.println("Analytical force = " + analiticalForce);
                System.out.println("Numerical force  = " + numericalForce);

                System.out.println("diff = " + diff + "\n" +
                        "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+EnergyElement.VERY_SMALL) = " +
                        2 * diff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + EnergyElement.VERY_SMALL));
                System.out.println();
            }
            if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\ne1 = " + e1);
            if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\ne2 = " + e2);
            if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\nanaliticalForce = " + analiticalForce);
        }

    }
}

