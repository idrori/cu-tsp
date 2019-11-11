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

public class TwoTorsionsCreator extends EnergyCreator implements KeyWords {
    protected ParametersList parametersList = null;

    public TwoTorsionsCreator(double weight) {
        super(weight,InfoType.TWO_TORSIONS);
        if (weight == 1)
            System.out.println("\n\n********************************************************" +
                    " TwoTorsion weight should not be one !! see TODO in TwoTorsionEnergyElement !" +
                    "******************************************************************************8\n\n");
    }

    public TwoTorsionsCreator() {
        super(InfoType.TWO_TORSIONS);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if (parametersList == null) {
            // Appanding the path to the filename list of the parameters.
            String[] strlist = new String[TWO_TORSIONS_PARAMETERS.length];
            String pathname = parametersDirectory(commands).concat("/");
            for (int cc = 0; cc < TWO_TORSIONS_PARAMETERS.length; cc++)
                strlist[cc] = pathname.concat(TWO_TORSIONS_PARAMETERS[cc]);

            parametersList = new TwoTorsionsParametersList(strlist);
        }
        TorsionPairList torsionPairList = TorsionPairList.createQuickAndDirtyTorsionPairList(protein, distanceMatrix);
        TorsionPairList relevantTorsionPairList = torsionPairList.filter(new HaveParametersFilter(parametersList));
        EnergyInfoElement info = new EnergyInfoElement(InfoType.TWO_TORSIONS, " Two torsions energy", weight);
        return new TwoTorsionsEnergy(relevantTorsionPairList, distanceMatrix, (TwoTorsionsParametersList) parametersList, info,
                "2Torsion");
    }

}

