/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.outOfPlane;

import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.info.InfoType;
import meshi.util.string.*;
import meshi.geometry.*;
import meshi.util.*;
import meshi.util.filters.*;


public class OutOfPlaneCreator extends SimpleEnergyCreator implements KeyWords {


    public OutOfPlaneCreator() {
        super(InfoType.OO_PLANE);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if (parametersList == null)
            parametersList = new OutOfPlaneParametersList(parametersDirectory(commands) +
                    "/" + OUT_OF_PLANE_PARAMETERS);
        TorsionList torsionList = TorsionList.createQuickAndDirtyTorsionList(protein, distanceMatrix);
        TorsionList OOPlist = (TorsionList) torsionList.filter(new TorsionList.FilterOOP());
        TorsionList relevantTorsionList = (TorsionList) OOPlist.filter(new HaveParametersFilter(parametersList));
        EnergyInfoElement info = new EnergyInfoElement(InfoType.OO_PLANE, " Out of Plane energy", weight);
        return new OutOfPlaneEnergy(distanceMatrix, relevantTorsionList, (OutOfPlaneParametersList) parametersList, info);
    }
}
