/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomCore;
import meshi.parameters.AtomType;


/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Feb 15, 2010
 * Time: 10:46:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZSummaPolarElement {
    private static final double Z_ALPHA = 10, Z_WEIGHT = 0.0000000001;
    protected final static int NumOfClusters = CooperativeZ5typesParameters.NumOfClusters;
    protected final static double[] expOfZScoredAtomEnergies = new double[NumOfClusters];
    protected final static double[] zScoredAtomEnergiesDerivatives = new double[NumOfClusters];
    protected final static double[] zScoredAtomEnergies = new double[NumOfClusters];
    protected final static double[] zPenalty = new double[NumOfClusters];
    protected final static double[] zPenaltyDerivative = new double[NumOfClusters];
    private double[] zEnergy, zEnergyDerivative;
    private double allButI;
    private int sizeOfMS;
    private final double z_weight;

    private static double zMinusZ_ALPHA;
    double summPlusEpsilon, summ, summDivSummPlusEpsilon;
    double[][] mean = new double[NumOfClusters][AtomType.values().length];
    double[][] std = new double[NumOfClusters][AtomType.values().length];
    double[][] proportion = new double[NumOfClusters][AtomType.values().length];
    double[] atomEnergies;
    int nAtoms;
    int initialPickNumber = 0;  //starts from 0 to 2
    double koef = 1;
    double EPSILON = 0.00001;

    static final double WIDTH = CooperativeSumma5PolarTypesEnergy.WIDTH;
    static final double HEIHT = CooperativeSumma5PolarTypesEnergy.HEIHT;
    static double M = (HEIHT + 1) * WIDTH;
    double weight;

    double[][] derivedExpOfZScoredAtomEnergies;
    double[] z = new double[NumOfClusters];
    double[] energyForEachND = new double[NumOfClusters],
            derivativeForEachND = new double[NumOfClusters];
    double energy;
    String name;

    CooperativeZSummaPolarElement(double[][] mean, double[][] std, double[][] proportion, int nAtoms, double[] atomEnergies, double allWeight, DistanceMatrix distanceMatrix, String name, double myWeight) {
        this.mean = mean;
        this.std = std;
        this.proportion = proportion;
        this.nAtoms = nAtoms;
        sizeOfMS = distanceMatrix.molecularSystem.size();
        if (nAtoms == 0)
            throw new RuntimeException("Zero number of atoms in " + name + this + ". Weird polarElement definition.");
        this.atomEnergies = atomEnergies;
        this.weight = allWeight * myWeight;
        z_weight = Z_WEIGHT * myWeight;
        this.name = name;
        derivedExpOfZScoredAtomEnergies = new double[NumOfClusters][sizeOfMS];
        zEnergy = new double[sizeOfMS];
        zEnergyDerivative = new double[sizeOfMS];


    }

    void resetZ() {
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            z[i] = 0;
            energyForEachND[i] = 0;
            derivativeForEachND[i] = 0;
            for (int k = 0; k < derivedExpOfZScoredAtomEnergies[i].length; k++)
                derivedExpOfZScoredAtomEnergies[i][k] = 0;
        }
        energy = 0;
    }

    void addZScore(AtomCore atom) {
        summ = 0;
        int number = atom.number;
        int type = atom.type().ordinal();
        zEnergy[number] = z_weight;
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            if (std[i][type] == 0)
                throw new RuntimeException("Zero std value for " + atom);
            zScoredAtomEnergies[i] = zScore(atomEnergies[number], mean[i][type], std[i][type]);
            if (zScoredAtomEnergies[i] > Z_ALPHA) {
                zMinusZ_ALPHA = zScoredAtomEnergies[i] - Z_ALPHA;
            } else if (zScoredAtomEnergies[i] < -Z_ALPHA) {
                zMinusZ_ALPHA = zScoredAtomEnergies[i] + Z_ALPHA;
            } else zMinusZ_ALPHA = 0;
            zPenalty[i] = zMinusZ_ALPHA * zMinusZ_ALPHA * zMinusZ_ALPHA * zMinusZ_ALPHA;
            zPenaltyDerivative[i] = 4 * zMinusZ_ALPHA * zMinusZ_ALPHA * zMinusZ_ALPHA / std[i][type];
            zEnergy[number] *= zPenalty[i];

            expOfZScoredAtomEnergies[i] = proportion[i][type] * Math.exp(-0.5 * koef * (zScoredAtomEnergies[i] * zScoredAtomEnergies[i]));

            summ += expOfZScoredAtomEnergies[i];
        }

//proceed Z[i] and derivatives
        summPlusEpsilon = summ + EPSILON;
        if (summPlusEpsilon == 0) throw new RuntimeException("weird summPlusEpsilon " + name + " " + summPlusEpsilon);
        summDivSummPlusEpsilon = summ / summPlusEpsilon;
        zEnergyDerivative[number] = 0;
        if (zEnergy[number] > 0) {
            for (int i = initialPickNumber; i < NumOfClusters; i++) {
                allButI = 1;
                for (int j = initialPickNumber; j < NumOfClusters; j++) {
                    if (i != j) allButI *= zPenalty[j];
                }
                zEnergyDerivative[number] += z_weight * zPenaltyDerivative[i] * allButI;
            }
        }


        double expEnergy, summNotI;
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            expEnergy = expOfZScoredAtomEnergies[i] / summPlusEpsilon;
            z[i] += (zScoredAtomEnergies[i] * expEnergy);
            summNotI = 0;
            for (int j = initialPickNumber; j < NumOfClusters; j++) {
                if (j != i)
                    summNotI += expOfZScoredAtomEnergies[j] * koef * (zScoredAtomEnergies[j] / std[j][type] - zScoredAtomEnergies[i] / std[i][type]);
            }
            summNotI += (-(1) * koef * zScoredAtomEnergies[i] / std[i][type] * EPSILON);
            derivedExpOfZScoredAtomEnergies[i][number] = expEnergy / std[i][type] + zScoredAtomEnergies[i] * (expOfZScoredAtomEnergies[i] * summNotI / (summPlusEpsilon * summPlusEpsilon));
        }

    }


    double zScore(double curE, double mean, double std) {
        if (std == 0)
            //               return 0;
            throw new RuntimeException("Zero std value in " + this);
        return (curE - mean) / std;
    }

    void makeEnergy() {
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            if ((nAtoms == 0) && (z[i] != 0))
                throw new RuntimeException("Something weird with number of atoms in " + name + this);
            if ((nAtoms != 0) && (z[i] == 0))
                throw new RuntimeException("Something weird with Z score calculation in " + name + this);
            z[i] /= nAtoms;
            calcEnergy2(i);
            //        calcEnergy4(i);
            //calcEnergy4withZI(i);
            if ((!(energyForEachND[i] < 0)) & (!(energyForEachND[i] == 0)) & (!(energyForEachND[i] > 0)))
                throw new RuntimeException("weird energy " + name + " " + energyForEachND[i]);
            energy += energyForEachND[i];
        }

        for (int i = 0; i < sizeOfMS; i++) {
            energy += zEnergy[i];
        }



    }

    void calcEnergy2(int i) {
        /* ------z^2 */
        double zz = z[i];
        double z2 = zz * zz;
        if (WIDTH != 10)
            throw new RuntimeException("Check the WIDTH value in CooperativeSumma5PolarTypesEnergy. It does not fit for this case in " + name + this);
        energyForEachND[i] = weight * (M / WIDTH - (M + z2) / (z2 + WIDTH));
        derivativeForEachND[i] = weight * (-2 * zz * (WIDTH - M) / (z2 + WIDTH) / (z2 + WIDTH) / nAtoms);
    }


    void evaluateAtoms(AtomCore atom1, AtomCore atom2, SummaAttribute summaAttribute, int order) {
        double derivative, force;
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            derivative = derivedExpOfZScoredAtomEnergies[i][atom2.number];
            //       if (derivative == 0.0)
            //            // throw new RuntimeException("Something is weird in Parameters "+derivative+name);
            //    System.out.println("Something is weird in Parameters "+derivative+name);

            force = derivativeForEachND[i] * derivative / 2;
            force = force * order;
            if (!atom1.status().frozen())
                atom1.addForce(force * summaAttribute.fx, force * summaAttribute.fy, force * summaAttribute.fz);
            if (!atom2.status().frozen())
                atom2.addForce(-force * summaAttribute.fx, -force * summaAttribute.fy, -force * summaAttribute.fz);
        }
        force = zEnergyDerivative[atom2.number] * order / 2.0;

        if ((!(force < 0)) & (!(force == 0)) & (!(force > 0)))
            throw new RuntimeException("this is weird1 " + force + " " + atom2);
        if (!atom1.status().frozen())
            atom1.addForce(force * summaAttribute.fx, force * summaAttribute.fy, force * summaAttribute.fz);
        if (!atom2.status().frozen())
            atom2.addForce(-force * summaAttribute.fx, -force * summaAttribute.fy, -force * summaAttribute.fz);

    }

}
