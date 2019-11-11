/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 05/08/2009
 * Time: 15:52:40
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZStdRamachandranCreator extends EnergyCreator {

    private static CooperativeZRamachandranParameters parameters = null;
    private RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator;
    private CooperativeZRamachandranCreator cooperativeZRamachandranCreator;

    public CooperativeZStdRamachandranCreator(CooperativeZRamachandranCreator cooperativeZRamachandranCreator) {
        super(InfoType.COOP_Z_STD_RAMACHANDRAN);
        this.cooperativeZRamachandranCreator = cooperativeZRamachandranCreator;
        this.ramachandranSidechainEnergyCreator = cooperativeZRamachandranCreator.ramachandranSidechainEnergyCreator;
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        /* load parameters if this is the first time the function is called */
        parameters = cooperativeZRamachandranCreator.parameters;
        if (parameters == null)
            throw new RuntimeException("No parameters for " + this);


//        RamachandranSidechainEnergy ramachandranSidechainEnergy = (RamachandranSidechainEnergy) ramachandranSidechainEnergyCreator.term();
        if (term != null) return term;
//        if (ramachandranSidechainEnergy == null)
//            throw new RuntimeException("Apparently an attempt to insantiate CooperativeZStdRamachandran " +
//                    "before RamachandranSidechainEnergy");
        return term = new CooperativeZStdRamachandran(ramachandranSidechainEnergyCreator,
                new EnergyInfoElement(infoType, "Cooperative energy term that enforce proper distribution of Ramachamdran energies", weight),
                parameters);
    }

}


