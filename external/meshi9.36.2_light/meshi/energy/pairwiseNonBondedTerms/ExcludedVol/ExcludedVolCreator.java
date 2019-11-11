/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;
import meshi.util.info.InfoType;

public class ExcludedVolCreator extends EnergyCreator implements KeyWords {


    private int type = 0; // allows us to run one of the EV scenario
    private double Rfac = 1.0;
    private Filter filter = null;
    private static ExcludedVolParametersList parametersList = null;

//    public ExcludedVolCreator(double weight, int type) {
//        super(weight);
//        this.type = type;
//    }
//
//    public ExcludedVolCreator(double weight) {
//        super(weight);
//    }

    public ExcludedVolCreator() {
        super(InfoType.EXCLUDED_VOL);
    }
    public ExcludedVolCreator(InfoType type) {
        super(type);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if (parametersList == null)
            parametersList = new ExcludedVolParametersList(parametersDirectory(commands) +
                    "/" + EXCLUDED_VOL_PARAMETERS);
        EnergyInfoElement info = new EnergyInfoElement(InfoType.EXCLUDED_VOL, "Excluded volume energy", weight);
        term = new ExcludedVol(distanceMatrix, parametersList, type, info, Rfac, filter);
        return term;
    }
}
