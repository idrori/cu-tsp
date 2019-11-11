/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity2DwithPP;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;

public class CompositePropensity2DEnergy
        extends SimpleEnergyTerm
        implements CompositeTorsionsDefinitions {

    public CompositePropensity2DEnergy() {
    }

    public CompositePropensity2DEnergy(
            ResidueTorsionsList residueTorsionsList,
            DistanceMatrix distanceMatrix,
            CompositePropensity2DParametersList cppl,
            EnergyInfoElement info,
            String comment) {
        super(toArray(distanceMatrix, residueTorsionsList), cppl, info);

        this.comment = comment;
        createElementsList(residueTorsionsList);
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        ResidueTorsions resTorsions =
                (ResidueTorsions) baseElement;
        CompositePropensity2DParameters cpp =
                (CompositePropensity2DParameters) parameters;

        return new CompositePropensity2DEnergyElement(
                resTorsions, cpp, weight);
    }

}
