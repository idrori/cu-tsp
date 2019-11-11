package meshi.energy.solvation;

import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomList;

/**
 * Created by chen on 15/06/2015.
 */
public class SolvationEnergySS extends SolvationEnergy {
    public SolvationEnergySS(AtomList atomList,
                           DistanceMatrix dm,
                           SolvateParametersList parameters,
                           EnergyInfoElement info) {
        super(atomList, dm, parameters,info);
    }
    public SolvationDistanceAttributeCreator solvationDistanceAttributeCreator() {
        return solvationDistanceAttributeCreatorSS;
    }
}
