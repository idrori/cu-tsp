/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.angle;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.info.InfoType;
import meshi.util.string.*;
import meshi.util.*;

public class AngleCreator extends SimpleEnergyCreator implements KeyWords {
    public AngleCreator() {
        super(InfoType.ANGLE);
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        AtomPairList bondList = (AtomPairList) protein.bonds().filter(new GoodBonds());
        AngleList angleList = new QuickAndDirtyAngleList(bondList, distanceMatrix);
        if (parametersList == null)
            parametersList = new AngleParametersList(parametersDirectory(commands) +
                    "/" + ANGLE_PARAMETERS);
        EnergyInfoElement info = new EnergyInfoElement(InfoType.ANGLE, "Angle energy", weight);
        return term = new AngleEnergy(angleList, distanceMatrix, (AngleParametersList) parametersList, info);
    }
}
