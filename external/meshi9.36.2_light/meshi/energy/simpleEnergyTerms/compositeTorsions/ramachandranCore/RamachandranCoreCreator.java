package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

/**
 * Created with IntelliJ IDEA.
 * User: tetyanam
 * Date: 05/02/14
 * Time: 12:45
 * To change this template use File | Settings | File Templates.
 */
public class RamachandranCoreCreator extends EnergyCreator implements KeyWords {
    protected RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator;
    protected RamachandranSidechainEnergy ramachandranSidechainEnergy;

    public RamachandranCoreCreator(RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator) {
        super(InfoType.RAMACHANDRAN_CORE);
        this.ramachandranSidechainEnergyCreator = ramachandranSidechainEnergyCreator;
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        boolean ssFlag;
        if (term != null) return term;
        /* load parameters if this is the first time the function is called */
        if (ramachandranSidechainEnergy == null)
                ramachandranSidechainEnergy = (RamachandranSidechainEnergy) ramachandranSidechainEnergyCreator.term();
        if (ramachandranSidechainEnergy == null)
            throw new RuntimeException("Apparently an attempt to create RamachandranCore " +
                    "before RamachandranSidechainEnergy");
        ssFlag = commands.firstWordFilter(infoType).secondWord("ss").thirdWord().equals("ON");
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.RAMACHANDRAN_CORE, "the term affects only on phi,psi in the ramachandran core", weight);
        return term = new RamachandranCoreEnergy(ramachandranSidechainEnergy, energyInfoElement,ssFlag,weight);
    }

/*
    public void updateEnergyTerm(AbstractEnergy ramachandranSidechainEnergy) {
        */
/* load parameters if this is the first time the function is called *//*

        if (ramachandranSidechainEnergy == null)
            throw new RuntimeException("Apparently an attempt to create RamachandranCore " +
                    "before RamachandranSidechainEnergy");
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.RAMACHANDRAN_CORE, "the term affects only on phi,psi in the ramachandran core", weight);
        ssFlag = commands.firstWordFilter(infoType).secondWord("ss")!= null;
        this.term = new RamachandranCoreEnergy((RamachandranSidechainEnergy)ramachandranSidechainEnergy, energyInfoElement,ssFlag, weight);
    }
*/

    public void updateRamachEnergyTerm(AbstractEnergy ramachandranSidechainEnergy) {
        this.ramachandranSidechainEnergy = (RamachandranSidechainEnergy) ramachandranSidechainEnergy;
    }

}
