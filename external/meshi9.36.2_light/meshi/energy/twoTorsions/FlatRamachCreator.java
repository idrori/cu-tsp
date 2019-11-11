/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.twoTorsions;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.ParametersList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionPairList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

public class FlatRamachCreator extends EnergyCreator implements KeyWords {
    protected ParametersList parametersList = null;

    public FlatRamachCreator(double weight) {
        super(weight,InfoType.FLAT_RAMACH);
    }

    public FlatRamachCreator() {
        super(InfoType.FLAT_RAMACH);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if (parametersList == null) {
            // Appanding the path to the filename list of the parameters.
            String[] strlist = new String[FLAT_RAMACH_PARAMETERS.length];
            String pathname = parametersDirectory(commands).concat("/");
            for (int cc = 0; cc < FLAT_RAMACH_PARAMETERS.length; cc++)
                strlist[cc] = pathname.concat(FLAT_RAMACH_PARAMETERS[cc]);

            parametersList = new TwoTorsionsParametersList(strlist);
        }

        TorsionPairList torsionPairList = TorsionPairList.createQuickAndDirtyTorsionPairList(protein, distanceMatrix);
        TorsionPairList relevantTorsionPairList = torsionPairList.filter(new HaveParametersFilter(parametersList));
        EnergyInfoElement info = new EnergyInfoElement(InfoType.FLAT_RAMACH, "Flat Ramachandran energy", weight);
        return term = new TwoTorsionsEnergy(relevantTorsionPairList, distanceMatrix, (TwoTorsionsParametersList) parametersList, info,
                "FlatRamachandran");
    }
}
