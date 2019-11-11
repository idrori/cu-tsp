package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * Created by chen on 03/04/2016.
 */
public class SwellCreator extends ExcludedVolCreator {
    public SwellCreator() {
        super(InfoType.SWELL);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        EnergyInfoElement info = new EnergyInfoElement(InfoType.SWELL, "Swell", weight);
        term = new Swell(distanceMatrix, info);
        return term;
    }

}
