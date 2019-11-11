/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

public class AtomicPairwisePMFSummaCreator extends
        EnergyCreator implements KeyWords {

    private DistanceMatrix distanceMatrix = null;
    private AtomicPairwisePMFSummaParameters parameters = null;

    public AtomicPairwisePMFSummaCreator() {
        super(InfoType.ATOMIC_PAIRWISE_PMF_SUMMA);
    }

//    public AtomicPairwisePMFSummaCreator(double weight) {
//        super(weight);
//    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if ((term != null) && (this.distanceMatrix == distanceMatrix)) return term; //chen3/3/12
        this.distanceMatrix = distanceMatrix;
        /* load parameters if this is the first time the function is called */
        if (parameters == null)
            parameters = new AtomicPairwisePMFSummaParameters(commands);
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.ATOMIC_PAIRWISE_PMF_SUMMA, " Pairwise energy term by Summa and Levitt", weight);
        energyInfoElement.getChildren().add(new EnergyInfoElement(InfoType.SUMMA_STD,"The distribution of energy over the model"));
        return term = new AtomicPairwisePMFSumma(distanceMatrix, energyInfoElement, parameters);
    }

}
