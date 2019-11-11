/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity2DwithPP;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.*;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

/**
 * The propensity 2D energy calculates the proensity of each amino acid to be in a specific position
 * at the (phi,psi) torsion space.
 */
public class CompositePropensity2DCreator
        extends EnergyCreator
        implements KeyWords, MeshiPotential {
    private static CompositePropensity2DParametersList parametersList = null;

    public CompositePropensity2DCreator(double weight) {
        super(weight,InfoType.PROPENSITY);
    }

    public CompositePropensity2DCreator() {
        super(1.0,InfoType.PROPENSITY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands) {
        /* retrieve parameters */
        if (parametersList == null) {
            String cpplFileName = parametersDirectory(commands) + "/" + COMPOSITE_PROPENSITY_2D_WITH_PP_PARAMETERS;
            parametersList = new CompositePropensity2DParametersList(cpplFileName);
        }

        /* create residue torsions list for protein */
        ResidueTorsionsList rtl = (new ResidueTorsionsList(protein, distanceMatrix)).filterPhiPsiResidues();

        /* return energy */
        EnergyInfoElement info = new EnergyInfoElement(InfoType.PROPENSITY, " I do not know what this energy term does", weight);
        return new CompositePropensity2DEnergy(rtl, distanceMatrix,
                (CompositePropensity2DParametersList) parametersList, info, "prop2D");
    }
}
