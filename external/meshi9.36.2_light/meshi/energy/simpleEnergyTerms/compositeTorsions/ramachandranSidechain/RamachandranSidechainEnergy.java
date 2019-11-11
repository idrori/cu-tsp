/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.parameters.ResidueType;
import meshi.util.Stat;
import meshi.util.Utils;

/**
 * A Ramachandran plot and sidechain torsions optimization energy term.
 * Statistical analysis of residue backbone and sidechain torsions in a
 * large database of residue data has been smoothed using
 * polynomial spline interpolation.
 * For a given residue the energy value approximates the percentage
 * of finding its current backbone and sidechain torsion angles.
 *
 * @author El-ad David Amir
 */
public class RamachandranSidechainEnergy
        extends SimpleEnergyTerm
        implements CompositeTorsionsDefinitions {

    protected static double[] sumPerResidueType = new double[ResidueType.values().length];
    protected static double[] sm2PerResidueType = new double[ResidueType.values().length];
    protected static int[] numberOfResiduesPerType = new int[ResidueType.values().length];
    private Stat stat;
    public static double[] sumPerResidueType() {
        return sumPerResidueType;
    }

    public static double[] sm2PerResidueType() {
        return sm2PerResidueType;
    }

    public static int[] numberOfResiduesPerType() {
        return numberOfResiduesPerType;
    }

    public static double[] residueEnergies;
    protected ResidueTorsionsList residueTorsionsList;

    public ResidueTorsionsList residueTorsionsList() {
        return residueTorsionsList;
    }

    public RamachandranSidechainEnergy() {
    }

    public RamachandranSidechainEnergy(
            ResidueTorsionsList residueTorsionsList,
            DistanceMatrix distanceMatrix,
            RamachandranSidechainParametersList rspl,
            EnergyInfoElement info,
            String comment) {
        super(toArray(distanceMatrix, residueTorsionsList), rspl, info);

        this.comment = comment;
        createElementsList(residueTorsionsList);
        for (int i = 0; i < ResidueType.values().length; i++)
            numberOfResiduesPerType[i] = 0;
        for (ResidueTorsions rt : residueTorsionsList)
            numberOfResiduesPerType[rt.getResidueType().ordinal()]++;
        this.residueTorsionsList = residueTorsionsList;
        stat = new Stat();
    }

    public RamachandranSidechainEnergyElement getHighestEnergyElement() throws EvaluationException{
        RamachandranSidechainEnergyElement out = null;
        double highestEnergy = -100000;
        for (EnergyElement element : elementsList) {
            double e = element.evaluate();
            if (e > highestEnergy) {
                out = (RamachandranSidechainEnergyElement) element;
                highestEnergy = e;
            }
        }
        return out;
    }

    public void scaleWeight(double factor) {
        weight *= factor;
        for (EnergyElement element : elementsList) {
            ((RamachandranSidechainEnergyElement) element).scaleWeight(factor);
        }
        Utils.println("Ramachandran side chain weight = " + weight);
    }

    public EnergyInfoElement evaluate() {
        stat.reset();
        for (int i = 0; i < ResidueType.values().length; i++) {
            sumPerResidueType[i] = 0;
            sm2PerResidueType[i] = 0;
        }
        for (ResidueTorsions r : residueTorsionsList)
            r.resetEnergy();

        double e, energy = 0;
        if (elementsList == null)
            throw new RuntimeException("this is very weird");
        for (EnergyElement energyElement : elementsList) {
            e = energyElement.evaluate();
            energy += e;
            stat.add(e);
        }
        if (info == null) throw new RuntimeException("null info in " + this);
        info.getChildren().get(0).setValue(stat.getStd());
        info.setValue(energy);
       return info;
    }

    public void evaluateAtoms() throws EvaluationException{
        for (int i = 0; i < ResidueType.values().length; i++) {
            sumPerResidueType[i] = 0;
            sm2PerResidueType[i] = 0;
        }
        for (ResidueTorsions r : residueTorsionsList)
            r.resetEnergy();

        super.evaluateAtoms();
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        ResidueTorsions resTorsions =
                (ResidueTorsions) baseElement;
        RamachandranSidechainParameters rsp =
                (RamachandranSidechainParameters) parameters;

        switch (NUM_SIDECHAIN_TORSIONS[resTorsions.getResidueType().ordinal()]) {
            case 0:
                return new RamachandranSidechainEnergyElementChi0(
                        resTorsions, rsp, weight);
            case 1:
                return new RamachandranSidechainEnergyElementChi1(
                        resTorsions, rsp, weight);
            case 2:
                return new RamachandranSidechainEnergyElementChi2(
                        resTorsions, rsp, weight);
            case 3:
                return new RamachandranSidechainEnergyElementChi3(
                        resTorsions, rsp, weight);
            case 4:
                return new RamachandranSidechainEnergyElementChi4(
                        resTorsions, rsp, weight);
            default:
                throw new RuntimeException("unidentified residue type");
        }
    }



}
