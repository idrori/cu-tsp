/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.geometry.DistanceMatrix;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfoElement;

public class RamachandranEnergy
        extends SimpleEnergyTerm
        implements CompositeTorsionsDefinitions {

    private ResidueTorsionsList residueTorsionsList;

    public RamachandranEnergy() {
    }

    public RamachandranEnergy(
            ResidueTorsionsList residueTorsionsList,
            DistanceMatrix distanceMatrix,
            RamachandranParametersList cppl,
            EnergyInfoElement info,
            String comment) {
        super(toArray(distanceMatrix, residueTorsionsList), cppl, info);

        this.comment = comment;
        createElementsList(residueTorsionsList);
        this.residueTorsionsList = residueTorsionsList;
    }

    public ResidueTorsionsList residueTorsionsList() {
        return residueTorsionsList;
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        ResidueTorsions resTorsions =
                (ResidueTorsions) baseElement;
        RamachandranParameters cpp =
                (RamachandranParameters) parameters;

        return new RamachandranEnergyElement(
                resTorsions, cpp, weight);
    }

    private void reset() {
        for (EnergyElement element : elementsList) {
            ((RamachandranEnergyElement) element).reset();
        }
    }

    public static void resetRamachandran(TotalEnergy totalEnergy) {
        RamachandranEnergy ramachandranEnergy = (RamachandranEnergy) totalEnergy.getEnergyTerm(new RamachandranEnergy());
        if (ramachandranEnergy != null) ramachandranEnergy.reset();
    }
}
