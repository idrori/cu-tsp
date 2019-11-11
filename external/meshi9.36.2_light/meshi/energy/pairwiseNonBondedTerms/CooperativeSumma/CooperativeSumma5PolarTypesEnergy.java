/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.energy.*;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.parameters.AtomType;
import meshi.parameters.ResidueType;
import meshi.util.MeshiAttribute;
import meshi.util.UpdateableException;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 16/02/2010
 * Time: 13:26:47
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeSumma5PolarTypesEnergy extends AbstractEnergy {

    private AtomicPairwisePMFSumma atomicPairwisePMFSumma;

    static final double WIDTH = 10; //z^2
    //static final double WIDTH = 30; //z^4
    //static final double WIDTH = 20;  //z^2 with zeroSpace
    static final double HEIHT = 1;

    private DistanceMatrix distanceMatrix;

    private SummaAttribute summaAttribute;
    private AtomCore atom1, atom2;
    private AtomType type1, type2;

    private CooperativeSummaProcesser cooperativeSummaProcesser;
    CooperativeZSummaPolarElement polarElement;
    CooperativeZSummaPolarElement nonPolarElement;
    CooperativeZSummaPolarElement neutralElement;
    CooperativeZSummaPolarElement polarBbElement;
    CooperativeZSummaPolarElement polarNN_OOElement;

//    public CooperativeSumma5PolarTypesEnergy() {}
public CooperativeSumma5PolarTypesEnergy(){}

    public CooperativeSumma5PolarTypesEnergy(DistanceMatrix distanceMatrix, EnergyInfoElement info, AtomicPairwisePMFSumma atomicPairwisePMFSumma,
                                             CooperativeZSummaPolarElement polarElement,
                                             CooperativeZSummaPolarElement nonPolarElement,
                                             CooperativeZSummaPolarElement neutralElement,
                                             CooperativeZSummaPolarElement polarNN_OOElement,
                                             CooperativeZSummaPolarElement polarBbElement,
                                             String comment) {
        super(toArray(distanceMatrix), info, EnergyType.NON_DIFFERENTIAL);
        if (weight != 0) {
            this.comment = comment;
            this.atomicPairwisePMFSumma = atomicPairwisePMFSumma;
            if (atomicPairwisePMFSumma == null)
                throw new RuntimeException("CooperativeAtomicPairwisePMFSumma term must not be created " +
                        "before the non-cooperative term is created.");
            cooperativeSummaProcesser = AtomicPairwisePMFSumma.cooperativeSummaProcesser;
            if (cooperativeSummaProcesser == null)
                throw new RuntimeException("The AtomicPairwisePMFSumma data have not been processed for " + this);
            this.distanceMatrix = distanceMatrix;

            this.polarElement = polarElement;
            this.nonPolarElement = nonPolarElement;
            this.neutralElement = neutralElement;
            this.polarBbElement = polarBbElement;
            this.polarNN_OOElement = polarNN_OOElement;
        }
    }


    public EnergyInfoElement evaluate() {
        if (on)
            return evaluate(false);
        else {
            info.setValue(0);
            return info;
        }
    }

    public void evaluateAtoms() {
        evaluate(true);
    }

    public EnergyInfoElement evaluate(boolean evaluateAtoms) {
        DistanceLists nonBondedLists = distanceMatrix.nonBondedList();
        CooperativeSummaInfoElement cooperativeSummaInfoElement = (CooperativeSummaInfoElement) info;
        double energy;
        polarElement.resetZ();
        nonPolarElement.resetZ();
        neutralElement.resetZ();
        if (polarBbElement != null) polarBbElement.resetZ();
        polarNN_OOElement.resetZ();
        for (AtomCore atom : distanceMatrix.molecularSystem) {
            if (atom.type().isOther())
                continue;
            if (atom.atom.residue().type == ResidueType.CYS) continue;
            if (atom.type().isHydrogen())
                throw new RuntimeException("Something weird with the H atom polarity (it must be OTHER)" + this);

            polarElement.addZScore(atom);
            nonPolarElement.addZScore(atom);
            neutralElement.addZScore(atom);
            if (atom.type().isPolarBackbone()) {
                if (polarBbElement != null) polarBbElement.addZScore(atom);
                polarNN_OOElement.addZScore(atom);
            }
        }

        polarElement.makeEnergy();
        nonPolarElement.makeEnergy();
        neutralElement.makeEnergy();

        if (polarBbElement != null) polarBbElement.makeEnergy();
        polarNN_OOElement.makeEnergy();

        cooperativeSummaInfoElement.COOPERATIVE_SUMMA_POLAR.setValue(polarElement.energy);
        cooperativeSummaInfoElement.COOPERATIVE_SUMMA_NON_POLAR.setValue(nonPolarElement.energy);
        cooperativeSummaInfoElement.COOPERATIVE_SUMMA_NEUTRAL.setValue(neutralElement.energy);
        cooperativeSummaInfoElement.COOPERATIVE_SUMMA_POLAR_NN_OO.setValue(polarNN_OOElement.energy);


        if (polarBbElement != null) {
            energy = polarElement.energy + nonPolarElement.energy + neutralElement.energy + polarBbElement.energy + polarNN_OOElement.energy;
            cooperativeSummaInfoElement.COOPERATIVE_SUMMA_POLAR_BB.setValue(polarBbElement.energy);
        } else {
            energy = polarElement.energy + nonPolarElement.energy + neutralElement.energy + polarNN_OOElement.energy;
            cooperativeSummaInfoElement.COOPERATIVE_SUMMA_POLAR_BB.setValue(0);
        }
        /*if (evaluateAtoms) {
            Utils.println(polarElement.name+"  "+polarElement.energy);
            Utils.println(nonPolarElement.name+"  "+nonPolarElement.energy);
            Utils.println(neutralElement.name+"  "+neutralElement.energy);
            Utils.println(polarNN_OOElement.name+"  "+polarNN_OOElement.energy);
            if (polarBbElement!=null) Utils.println(polarBbElement.name+"  "+polarBbElement.energy);
            //for (AtomCore atom:distanceMatrix.molecularSystem){
              //  atom.atom.addEnergy(energy/distanceMatrix.molecularSystem.size());
              // }
        }*/


        for (DistanceList distances : nonBondedLists) {
            for (Distance distance :distances) {
                if (distance.mode().dead) continue;
                if (distance.mode().frozen) continue;
                summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
                if (summaAttribute == null) continue;
                atom1 = distance.atom1;
                atom2 = distance.atom2;
                type1 = atom1.type();
                type2 = atom2.type();


                if (type1.isOther() || type2.isOther())
                    continue;

                if (type1.isPolarBackbone() && type2.isPolarBackbone()) {

                    if ((type1.isOxygen() && type2.isOxygen()) ||
                            (type1.isNitrogen() && type2.isNitrogen())) {
                        polarNN_OOElement.evaluateAtoms(atom1, atom2, summaAttribute, 1);
                        polarNN_OOElement.evaluateAtoms(atom2, atom1, summaAttribute, -1);
                    } else {
                        if (Math.abs(atom1.atom.residue().number() - atom2.atom.residue().number()) < 3)
                            continue;
                        if (polarBbElement != null) polarBbElement.evaluateAtoms(atom1, atom2, summaAttribute, 1);
                        if (polarBbElement != null) polarBbElement.evaluateAtoms(atom2, atom1, summaAttribute, -1);
                    }
                } else {
                    if (type1.isPolar())
                        polarElement.evaluateAtoms(atom1, atom2, summaAttribute, 1);
                    else if (type1.isNonPolar())
                        nonPolarElement.evaluateAtoms(atom1, atom2, summaAttribute, 1);
                    else if (type1.isNeutral())
                        neutralElement.evaluateAtoms(atom1, atom2, summaAttribute, 1);
                    else
                        throw new RuntimeException("Some weird type of the atom polarity (see AtomTypes.java) in " + this);

                    if (type2.isPolar())
                        polarElement.evaluateAtoms(atom2, atom1, summaAttribute, -1);
                    else if (type2.isNonPolar())
                        nonPolarElement.evaluateAtoms(atom2, atom1, summaAttribute, -1);
                    else if (type2.isNeutral())
                        neutralElement.evaluateAtoms(atom2, atom1, summaAttribute, -1);
                    else
                        throw new RuntimeException("Some weird type of the atom polarity (see AtomTypes.java) in " + this);
                }
            }
        }

        info.setValue(energy);
      /* if ((energy > 1E18) || (energy < -1E18))  {
             System.out.println(polarElement.energy + "     "+nonPolarElement.energy +"     " +
                                                     neutralElement.energy + "      "+ polarBbElement.energy + "       "+polarNN_OOElement.energy);
            throw new RuntimeException("Bad energy "+this+"  "+energy+"    "+info);
        }     */
          return info;
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
            atomicPairwisePMFSumma.evaluate();
            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            System.out.println(" evaluating e1 " + i);
            double e1 = evaluate(true).energy();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += EnergyElement.DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
                totalEnergy.update();
            System.out.println(" evaluating e2 " + i + " " + EnergyElement.DX);
            atomicPairwisePMFSumma.evaluate();
            double e2 = evaluate(true).energy();
            double de = e2 - e1;
            double numericalForce = -de / EnergyElement.DX;
            coordinates[i][0] -= EnergyElement.DX;
                totalEnergy.update();

            double diff = Math.abs(analiticalForce - numericalForce);

            if ((2 * diff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + EnergyElement.VERY_SMALL)) > EnergyElement.relativeDiffTolerance) {
                System.out.println("Testing " + this);
                System.out.println("Atom[" + atom + "]." + EnergyElement.XYZ.charAt(i) + " = " + x);
                System.out.println("Analytical force = " + analiticalForce);
                System.out.println("Numerical force  = " + numericalForce);

                System.out.println("diff = " + diff + "\n" +
                        "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+EnergyElement.VERY_SMALL) = " +
                        2 * diff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + EnergyElement.VERY_SMALL));
            }
            if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\ne1 = " + e1);
            if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\nz2 = " + e2);
            if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\nanaliticalForce = " + analiticalForce);
        }

    }
}


