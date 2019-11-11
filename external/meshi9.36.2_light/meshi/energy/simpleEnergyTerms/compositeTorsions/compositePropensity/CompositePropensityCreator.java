/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.*;
import meshi.util.*;
import meshi.util.filters.*;
import meshi.util.info.InfoType;

/**
 * The propensity energy identifies each amino acid's favorable positions
 * reset the (phi,psi,chi_1) torsion space. Alanine, Glycine and Proline do
 * not have a propensity value.
 */
public class CompositePropensityCreator
        extends EnergyCreator
        implements KeyWords, MeshiPotential {

    public CompositePropensityCreator() {
        super(InfoType.PROPENSITY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands) {
        /* retrieve parameters */
        String cpplFileName =
                parametersDirectory(commands) + "/" + COMPOSITE_PROPENSITY_PARAMETERS;
        CompositePropensityParametersList cppl =
                new CompositePropensityParametersList(cpplFileName);

        /* create residue torsions list for protein */
        ResidueTorsionsList rtl = new ResidueTorsionsList(
                protein, distanceMatrix);

        /* return energy */
        EnergyInfoElement info = new EnergyInfoElement(InfoType.PROPENSITY, "Torsion propensity energy", weight);
        return term = new CompositePropensityEnergy(
                rtl, distanceMatrix, cppl, info, "torsionProp");
    }
}
