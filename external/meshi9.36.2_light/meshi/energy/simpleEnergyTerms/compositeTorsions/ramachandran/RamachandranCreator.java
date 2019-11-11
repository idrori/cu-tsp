/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;
import meshi.parameters.*;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

/**
 * The Ramachandran energy calculates the probability of each amino acid to be in a specific position
 * at the (phi,psi) torsion space.
 */
public class RamachandranCreator extends EnergyCreator implements KeyWords, MeshiPotential {
    private static RamachandranParametersList parametersList = null;

    public RamachandranCreator(double weight) {
        super(weight,InfoType.RAMACHANDRAN);
    }

    public RamachandranCreator() {
        super(InfoType.RAMACHANDRAN);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands) {
        /* retrieve parameters */
        if (parametersList == null) {
            String cpplFileName = parametersDirectory(commands) + "/" + COMPOSITE_PROPENSITY_2D_PARAMETERS;
            parametersList = new RamachandranParametersList(cpplFileName);
        }

        /* create residue torsions list for protein */
        ResidueTorsionsList rtl = (new ResidueTorsionsList(protein, distanceMatrix)).filterPhiPsiResidues();

        /* return energy */
        EnergyInfoElement info = new EnergyInfoElement(InfoType.RAMACHANDRAN, "Ramachandran energy", weight);
        return term = new RamachandranEnergy(rtl, distanceMatrix,
                parametersList, info, "Ramach");
    }
}
