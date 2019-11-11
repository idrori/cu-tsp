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
 * Date: 09/06/2009
 * Time: 13:21:12
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZRamachandranCreator extends EnergyCreator implements KeyWords {

    protected static CooperativeZRamachandranParameters parameters = null;
    protected RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator;

    public CooperativeZRamachandranCreator(RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator) {
        super(InfoType.COOP_Z_RAMACHANDRAN);
        this.ramachandranSidechainEnergyCreator = ramachandranSidechainEnergyCreator;
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        if (term != null) return term;
        /* load parameters if this is the first time the function is called */
        if (parameters == null)
            parameters = new CooperativeZRamachandranParameters(commands);
//        RamachandranSidechainEnergy ramachandranSidechainEnergy = (RamachandranSidechainEnergy) ramachandranSidechainEnergyCreator.term();
//        if (ramachandranSidechainEnergy == null)
//            throw new RuntimeException("Apparently an attempt to insantiate CooperativeZRamachandran " +
//                    "before RamachandranSidechainEnergy");
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.COOP_Z_RAMACHANDRAN, "A cooperative term that affect the distribution of Ramachandran_sidchain energies", weight);
        return term = new CooperativeZRamachandran(ramachandranSidechainEnergyCreator, energyInfoElement, parameters);
    }

}

