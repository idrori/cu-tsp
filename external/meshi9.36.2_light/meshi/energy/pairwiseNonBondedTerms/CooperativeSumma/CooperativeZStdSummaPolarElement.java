/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.molecularElements.atoms.AtomCore;
import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 16/02/2010
 * Time: 10:58:20
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZStdSummaPolarElement extends CooperativeZSummaPolarElement {
    private double[] z2sum = new double[NumOfClusters];
    // double [] expOfZScoredAtomEnergies = new double[NumOfClusters];
    //double [] zScoredAtomEnergies = new double[NumOfClusters];

    private double[][] energyPerAtom;
    double stdProtein;
    double stdProteinMinus1;
    double stdProteinMinus1in2;


    CooperativeZStdSummaPolarElement(double[][] mean, double[][] std, double[][] proportion, int nAtoms, double[] atomEnergies, double allWeight, DistanceMatrix distanceMatrix, String name, double myWeight) {
        super(mean, std, proportion, nAtoms, atomEnergies, allWeight, distanceMatrix, name, myWeight);
        energyPerAtom = new double[NumOfClusters][distanceMatrix.molecularSystem.size()];
    }

    void resetZ() {
        super.resetZ();
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            z2sum[i] = 0;
            for (int k = 0; k < energyPerAtom[i].length; k++)
                energyPerAtom[i][k] = 0;
        }
    }

    void addZScore(AtomCore atom) {
        //if (atomEnergies[atom.number]==0) return;
        double summ = 0;
        int number = atom.number;
        int type = atom.type().ordinal();

        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            if (std[i][type] == 0)
                throw new RuntimeException("Zero std value for " + atom);
            //continue;
            zScoredAtomEnergies[i] = zScore(atomEnergies[number], mean[i][type], std[i][type]);
            expOfZScoredAtomEnergies[i] = proportion[i][type] * Math.exp(-0.5 * koef * (zScoredAtomEnergies[i] * zScoredAtomEnergies[i]));
            summ += expOfZScoredAtomEnergies[i];
        }
        summ += EPSILON;
//proceed Z[i] and derivatives
        if (summ == 0) //throw new RuntimeException("weird summPlusEpsilon "+name+" "+summPlusEpsilon);
            System.out.println("weird summPlusEpsilon in " + name + " " + summ);
        if (summ != 0) {
            double ze;
            double summNotI;
            double expEnergyPerAtom;
            for (int i = initialPickNumber; i < NumOfClusters; i++) {
                expEnergyPerAtom = expOfZScoredAtomEnergies[i] / summ;
                ze = zScoredAtomEnergies[i] * expEnergyPerAtom;
                energyPerAtom[i][number] = ze;
                z[i] += ze;
                z2sum[i] += (ze * ze);

                summNotI = 0;
                for (int j = initialPickNumber; j < NumOfClusters; j++) {
                    if (j != i)
                        summNotI += expOfZScoredAtomEnergies[j] * koef * (zScoredAtomEnergies[j] / std[j][type] - zScoredAtomEnergies[i] / std[i][type]);
                }
                summNotI += (-(1) * koef * zScoredAtomEnergies[i] / std[i][type] * EPSILON);
                derivedExpOfZScoredAtomEnergies[i][number] = expEnergyPerAtom / std[i][type] + zScoredAtomEnergies[i] * (expOfZScoredAtomEnergies[i] * summNotI / (summ * summ));
            }
        } else {
            for (int i = initialPickNumber; i < NumOfClusters; i++) {
                energyPerAtom[i][number] = 0;
                derivedExpOfZScoredAtomEnergies[i][number] = 0;
            }
        }

    }

    void calcEnergy2(int i) {
        /* ------stdProteinMinus1^2 */
        stdProtein = z2sum[i] / nAtoms - z[i] * z[i];
        stdProteinMinus1 = stdProtein - 1;
        stdProteinMinus1in2 = stdProteinMinus1 * stdProteinMinus1;
        if (WIDTH != 10)
            throw new RuntimeException("Check the WIDTH value in CooperativeSumma5PolarTypesEnergy. It does not fit for this case in " + this);
        energyForEachND[i] = weight * (M / WIDTH - (M + stdProteinMinus1in2) / (stdProteinMinus1in2 + WIDTH));
        derivativeForEachND[i] = weight * (4 * stdProteinMinus1 * (M - WIDTH) / (stdProteinMinus1in2 + WIDTH) / (stdProteinMinus1in2 + WIDTH) / nAtoms);
    }

    void calcEnergy4(int i) {
        /* ------stdProteinMinus1^4 */
        stdProtein = z2sum[i] / nAtoms - z[i] * z[i];
        stdProteinMinus1 = stdProtein - 1;
        stdProteinMinus1in2 = stdProteinMinus1 * stdProteinMinus1;
        double stdProteinMinus1in4 = stdProteinMinus1in2 * stdProteinMinus1in2;
        if (WIDTH != 30)
            throw new RuntimeException("Check the WIDTH value in CooperativeSumma5PolarTypesEnergy. It does not fit for this case in " + this);
        energyForEachND[i] = weight * (M / WIDTH - (M + stdProteinMinus1in4) / (stdProteinMinus1in4 + WIDTH));
        derivativeForEachND[i] = weight * 2 * 4 * stdProteinMinus1 * stdProteinMinus1in2 * (M - WIDTH) / (stdProteinMinus1in4 + WIDTH) / (stdProteinMinus1in4 + WIDTH) / nAtoms;
    }

    void calcEnergy4withZI(int i) {
        /* ------stdProteinMinus1^4  with Zero Interval*/
        stdProtein = z2sum[i] / nAtoms - z[i] * z[i];
        stdProteinMinus1 = stdProtein - 1;
        int coef;
        double zeroField = 0.15;
        double stdMinusZF;
        if ((stdProteinMinus1 >= -zeroField) && (stdProteinMinus1 <= zeroField)) {
            energyForEachND[i] = 0;
            derivativeForEachND[i] = 0;
        } else {
            if (stdProteinMinus1 < -zeroField)
                coef = -1;
            else coef = 1;                          // (z > zeroField)
            stdMinusZF = stdProteinMinus1 - coef * zeroField;
            stdProteinMinus1in2 = stdMinusZF * stdMinusZF;
            if (WIDTH != 20)
                throw new RuntimeException("Check the WIDTH value = " + WIDTH + " in CooperativeSumma5PolarTypesEnergy. It does not fit for this case in " + this);
            double stdProteinMinus1in4 = stdProteinMinus1in2 * stdProteinMinus1in2;
            energyForEachND[i] = weight * (M / WIDTH - (M + stdProteinMinus1in4) / (stdProteinMinus1in4 + WIDTH));
            derivativeForEachND[i] = weight * 2 * 4 * stdMinusZF * stdProteinMinus1in2 * (M - WIDTH) / (stdProteinMinus1in4 + WIDTH) / (stdProteinMinus1in4 + WIDTH) / nAtoms;
        }
    }

    void evaluateAtoms(AtomCore atom1, AtomCore atom2, SummaAttribute summaAttribute, int order) {
        double derivative, force;
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            derivative = derivedExpOfZScoredAtomEnergies[i][atom2.number];

            force = derivativeForEachND[i] * derivative * (energyPerAtom[i][atom2.number] - z[i]) / 2;
            force = force * order;
            if (!atom1.status().frozen())
                atom1.addForce(force * summaAttribute.fx, force * summaAttribute.fy, force * summaAttribute.fz);
            if (!atom2.status().frozen())
                atom2.addForce(-force * summaAttribute.fx, -force * summaAttribute.fy, -force * summaAttribute.fz);
        }
    }

}
