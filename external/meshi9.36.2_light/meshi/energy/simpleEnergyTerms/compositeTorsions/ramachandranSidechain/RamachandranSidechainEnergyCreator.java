/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

public class RamachandranSidechainEnergyCreator
        extends EnergyCreator
        implements KeyWords, MeshiPotential {

    public RamachandranSidechainEnergyCreator() {
        super(InfoType.RAMACHANDRAN_SIDECHAIN);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix distanceMatrix, CommandList commands) {
        if (term != null) {
            if (term.hasTheSameDistanceMatrix(distanceMatrix))
            return term;
        }
        /* retrieve parameters */
        String rsplFileName =
                parametersDirectory(commands) + "/" + COMPOSITE_TORSIONS_PARAMETERS;
        RamachandranSidechainParametersList rspl =
                new RamachandranSidechainParametersList(rsplFileName);

        /* create residue torsions list for protein */
        ResidueTorsionsList rtl = new ResidueTorsionsList(
                protein, distanceMatrix);

        /* return energy */
        EnergyInfoElement info = new EnergyInfoElement(InfoType.RAMACHANDRAN_SIDECHAIN, "Torsion energy that combines Ramachandran and rotamers ", weight);
        info.getChildren().add(new EnergyInfoElement(InfoType.RAMACH_STD,"The distribution of energy over the model"));
        return term = new RamachandranSidechainEnergy(rtl, distanceMatrix, rspl, info, "EnResidue");
    }

}
