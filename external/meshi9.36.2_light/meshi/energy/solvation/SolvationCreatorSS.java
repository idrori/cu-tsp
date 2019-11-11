package meshi.energy.solvation;

import meshi.energy.AbstractEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.info.InfoType;

/**
 * Created by chen on 15/06/2015.
 */
public class SolvationCreatorSS extends SolvationCreator {
    public SolvationCreatorSS(double weight) {
        super(weight);
    }

    protected AbstractEnergy createTerm(AtomList atomList, DistanceMatrix distanceMatrix, SolvateParametersList parametersList) {
        SolvationInfoElement info = new SolvationInfoElement(InfoType.SOLVATION, "Solvation energy",
                weightSCPolarSolvate,
                weightSCCarbonSolvate,
                weightBBPolarSolvate,
                weightBBCarbonSolvation,
                weightStdSolvation);
        return new SolvationEnergySS(atomList,
                distanceMatrix,
                parametersList,
                info);
    }

}
