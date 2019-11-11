/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

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

public class CompositePropensityEnergy extends SimpleEnergyTerm implements CompositeTorsionsDefinitions {
    protected static double[] sumPerResidueType = new double[ResidueType.values().length];
    protected static double[] sm2PerResidueType = new double[ResidueType.values().length];
    protected static int[] numberOfResiduesPerType = new int[ResidueType.values().length];

    public static double[] sumPerResidueType() {
        return sumPerResidueType;
    }

    public static double[] sm2PerResidueType() {
        return sm2PerResidueType;
    }

    public static int[] numberOfResiduesPerType() {
        return numberOfResiduesPerType;
    }

    protected ResidueTorsionsList residueTorsionsList;

    public ResidueTorsionsList residueTorsionsList() {
        return residueTorsionsList;
    }

    public CompositePropensityEnergy() {
    }

    public CompositePropensityEnergy(
            ResidueTorsionsList residueTorsionsList,
            DistanceMatrix distanceMatrix,
            CompositePropensityParametersList cppl,
            EnergyInfoElement info,
            String comment) {
        super(toArray(distanceMatrix, residueTorsionsList), cppl, info);
        this.residueTorsionsList = residueTorsionsList;
        for (int i = 0; i < ResidueType.values().length; i++)
            numberOfResiduesPerType[i] = 0;

        this.comment = comment;
        createElementsList(residueTorsionsList);
    }

    public EnergyInfoElement evaluate() {
        for (int i = 0; i < ResidueType.values().length; i++) {
            sumPerResidueType[i] = 0;
            sm2PerResidueType[i] = 0;
        }
        for (ResidueTorsions r : residueTorsionsList)
            r.resetEnergy();

        try {
            return super.evaluate();
        } catch (EvaluationException ex) {throw new RuntimeException(ex);}
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        ResidueTorsions resTorsions =
                (ResidueTorsions) baseElement;
        CompositePropensityParameters cpp =
                (CompositePropensityParameters) parameters;

        return new CompositePropensityEnergyElement(
                resTorsions, cpp, weight);
    }

}
