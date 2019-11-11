/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.bond;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.util.info.InfoType;
import meshi.util.string.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.*;

/**
 * A factory for BondEnergy objects.
 */
public class BondCreator extends SimpleEnergyCreator implements KeyWords {
    private boolean simpleFlag;
    public BondCreator() {
        super(InfoType.BOND);
        simpleFlag = false;
    }


    public BondCreator(String s) {
        this();
        if (!s.equals("simple"))
            throw new RuntimeException("\"simple\" is the only supported string here");
        simpleFlag = true;
    }
    /**
     * <pre>
     * hides all the hard work needed to generate a BondEnergy object.
     * a) Extract the bonds from the protein.
     * b) Finds and reads the parameters file.
     */
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        AtomPairList bondList = (AtomPairList) protein.bonds().filter(new GoodBonds());
        if (parametersList == null)
            parametersList = new BondParametersList(parametersDirectory(commands) +
                    "/" + BOND_PARAMETERS);
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.BOND, "Bond energy", weight);
        term = new BondEnergy(bondList, distanceMatrix, (BondParametersList) parametersList, energyInfoElement,simpleFlag);
        return term;
    }
}
