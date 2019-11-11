/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.beta;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 18/01/2010
 * Time: 13:51:05
 * To change this template use File | Settings | File Templates.
 */
public class BetaCreator extends EnergyCreator implements KeyWords {
    public BetaCreator() {
        super(InfoType.BETA);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.BETA, "Beta sheet promoting term", weight);
        term = new BetaEnergy(protein, energyInfoElement);
        return term;
    }
}
