/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.energy.*;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeSummaProcesser;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.SummaAttribute;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.AtomStatus;
import meshi.util.MeshiAttribute;
import meshi.util.ResidueData;
import meshi.util.Stat;
import meshi.util.UpdateableException;
import meshi.util.info.ChainsInfo;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfo;

import java.util.ArrayList;

/**
 * Atomic pairwise potential of mean force, original energy term by Chris
 * Summa.
 * <p/>
 * Reference:
 * Summa CM and Levitt M, Near-native structure refinement using in vacuo
 * energy minimization. PNAS 2007;104(9):3177-82.
 *
 * @author El-ad David Amir
 */
public class AtomicPairwisePMFSumma extends AbstractEnergy implements EvaluatesResidues{
    public final DistanceMatrix distanceMatrix;
    private boolean frozenFlag = false;
    private AtomCore atom1, atom2;
    private int bin;
    private double d, energy0, energy1, x, x2, x3, fx, fy, fz, halfEnergy0, e;
    private int numberOfClashes;
    private ArrayList<Distance> clashes = new ArrayList<Distance>();
    private MeshiInfo stdElement;
    private Stat stat;
    public int numberOfClashes() {
        return numberOfClashes;
    }
    public String clashes() {
        String out = "clashes:";
        for (Distance distance : clashes)
            out += distance+"\n";
        return out;
    }
    /* polynomial spline for this pair */
    private static CoefficientsMatrixForAtomPairSpline coefs;

    //  private static int    numberOfTypes = AtomType.values().length;
    private AtomicPairwisePMFSummaParameters parameters;
    //protected static final double[] sumPerAtomType      = new double[numberOfTypes];
    //protected static final double[] sm2PerAtomType      = new double[numberOfTypes];
    //protected static       int[]    numberOfEachAtomType = null;
    //public    static final double[] sumPerAtomType() { return sumPerAtomType;}
    //public    static final double[] sm2PerAtomType() { return sm2PerAtomType;}
    //public    static       int[]    numberOfEachAtomType() {return numberOfEachAtomType;}
    public static int atom1TypeOrdinal, atom2TypeOrdinal, atomTypeOrdinal;

    public static CooperativeSummaProcesser cooperativeSummaProcesser;

    public AtomicPairwisePMFSumma() {
        distanceMatrix = null;
    }

    public AtomicPairwisePMFSumma(DistanceMatrix distanceMatrix, EnergyInfoElement info, AtomicPairwisePMFSummaParameters parameters) {
        super(toArray(distanceMatrix), info);
        this.distanceMatrix = distanceMatrix;
        comment = "AtomicPairwisePMFSumma";

        for (AtomCore atom : distanceMatrix.molecularSystem)
            if (atom.status() == AtomStatus.FROZEN) frozenFlag = true;
        this.parameters = parameters;

        cooperativeSummaProcesser = new CooperativeSummaProcesser(distanceMatrix);
        stdElement = info.getChildren().get(0);
        stat = new Stat();
    }

    public boolean evaluatesResidues() {return true;}
    /**
     * Evaluates energy for each distance
     *
     * @return a sum of all energy elements
     */
    public EnergyInfoElement evaluate() {
        return evaluateResidues(false, null);
    }


    /**
     * Describe <code>evaluateAtoms</code> method here.
     */
    public void evaluateAtoms() {
        evaluateResidues(true, null);

    }

    public EnergyInfoElement evaluate(boolean evaluateAtoms) {
        evaluateResidues(evaluateAtoms, null);
        return info;
    }

    public void evaluateResidues(ChainsInfo chainsInfo) {
        evaluateResidues(false, chainsInfo);
    }

    public EnergyInfoElement evaluateResidues(boolean evaluateAtoms, ChainsInfo chainsInfo) {

            double e;
        numberOfClashes = 0;
        clashes.clear();
        stat.reset();

        double energy = 0;
        DistanceLists nonBondedLists = distanceMatrix.nonBondedList();
        if (nonBondedLists == null) throw new RuntimeException("nonBondedList == null");

        //if (cooperativeSummaProcesser != null)
        cooperativeSummaProcesser.reset();

        ResidueData residueData = null;
        if (chainsInfo != null)
            residueData = new ResidueData(chainsInfo);

        SummaAttribute summaAttribute;
        for (DistanceList distanceList : nonBondedLists) {
            for (Distance distance : distanceList) {
                if (!distance.mode().frozen) {
                    distance.summaValue = 0;
                    //if (cooperativeSummaProcesser != null) {
                    summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
                    if (summaAttribute != null)
                        summaAttribute.reset();
                    //   }
                    e = evaluateElement(distance, evaluateAtoms);
                    if (chainsInfo != null)  {
                        residueData.add(distance.atom1(), e/2);
                        residueData.add(distance.atom2(), e/2);
                    }
                    if ((e > 8) && (Math.abs(distance.atom1().residueNumber() - distance.atom2().residueNumber()) > 1)) {
                        numberOfClashes++;
                        clashes.add(distance);
                    }
                    energy += e;
                    stat.add(e);
                }
            }
        }
        if (chainsInfo != null)
            chainsInfo.add(residueData, info.type);
        info.setValue(energy);
        stdElement.setValue(stat.getStd());
        return info;
    }


    public double evaluateElement(Distance distance, boolean evaluateAtoms) {
        atom1 = distance.atom1;
        atom2 = distance.atom2;
        atom1TypeOrdinal = atom1.type().ordinal();
        atom2TypeOrdinal = atom2.type().ordinal();
        coefs = parameters.splineForAtomPair(atom1TypeOrdinal, atom2TypeOrdinal);
        /* turn off term if any of the atoms was not found in parameters file */
        if (coefs == null) return 0;

        /* calculate energy and derivative */
        d = distance.distance();
        bin = (int) (d * 10);
// 		bin = coefs.bins.length-1;
// 		while( (bin > 0) && coefs.bins[bin] >= d )
// 			bin--;

        /* correct x according to its distance from start of bin */
        x = d - coefs.bins[bin];
        /* calculate powers of x */
        x2 = x * x;
        x3 = x2 * x;

        /* calculate energy */
        try {
            energy0 = coefs.coefs[bin][0] * x3 + coefs.coefs[bin][1] * x2 + coefs.coefs[bin][2] * x + coefs.coefs[bin][3];
        }
        catch (RuntimeException ex) {
            System.out.println("Problem in evaluate\n" +
                    atom1 + " " + atom1.type() + "\n" +
                    atom2 + " " + atom2.type() + "\n" +
                    "bin = " + bin + "\n" +
                    "distance = " + d);
            System.out.println("coefs = " + coefs);
            System.out.println("coefs.coefs = " + coefs.coefs);
            System.out.println("coefs.coefs = " + coefs.coefs);
            throw ex;
        }
        /* calculate first derivative */
        energy1 = 3 * coefs.coefs[bin][0] * x2 + 2 * coefs.coefs[bin][1] * x + coefs.coefs[bin][2];

        halfEnergy0 = energy0 * 0.5;
        distance.summaValue += energy0;

/* apply weight */
        energy0 = energy0 * weight;
        //energy1     = energy1 * weight;
        fx = -energy1 * distance.dx * distance.invDistance;
        fy = -energy1 * distance.dy * distance.invDistance;
        fz = -energy1 * distance.dz * distance.invDistance;
        //     if (cooperativeSummaProcesser != null)
        cooperativeSummaProcesser.setCooperativeStatistic(distance, fx, fy, fz, halfEnergy0);

        fx *= weight;
        fy *= weight;
        fz *= weight;

        if (frozenFlag) {
            if (!atom1.status().frozen()) atom1.addForce(fx, fy, fz);
            if (!atom2.status().frozen()) atom2.addForce(-fx, -fy, -fz);
        } else {
            atom1.addForce(fx, fy, fz);
            atom2.addForce(-fx, -fy, -fz);
        }

        if (evaluateAtoms) {
            atom1.atom.addEnergy(halfEnergy0);
            atom2.atom.addEnergy(halfEnergy0);
        }


        return energy0;
    }


    public void test(TotalEnergy totalEnergy, Atom atom) {
//	 System.out.println("Sorry no test for  AtomicPairwisePMFSumma yet");
        if (!on) {
            System.out.println("" + this + " is off");
            return;
        }
        System.out.println("Testing " + this + " using " + atom);
        if (atom == null)
            throw new RuntimeException("Cannot test " + this);

        DistanceLists nonBondedLists = distanceMatrix.nonBondedList();
        if (nonBondedLists == null) throw new RuntimeException("nonBondedList == null");
        for (DistanceList distances : nonBondedLists) {
            for (Distance distance : distances) {
                if (distance.mode().frozen) continue;
                if ((distance.atom1() != atom) && (distance.atom2() != atom)) continue;

                //energyElement test
                double[][] coordinates = new double[3][];
                coordinates[0] = atom.X();
                coordinates[1] = atom.Y();
                coordinates[2] = atom.Z();
                for (int i = 0; i < 3; i++) {
                        totalEnergy.update();
                    double x = coordinates[i][0];
                    coordinates[i][1] = 0;
                    double e1 = evaluateElement(distance, false);
                    double analiticalForce = coordinates[i][1];
                    coordinates[i][0] += EnergyElement.DX;
                    // Whatever should be updated ( such as distance matrix torsion list etc. )
                        totalEnergy.update();
                    double e2 = evaluateElement(distance, false);
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
    }

}
