/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 09/06/2009
 * Time: 10:08:17
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZPropensityCreator extends EnergyCreator implements KeyWords {

    protected static CooperativeZPropensityParameters parameters = null;
    protected CompositePropensityCreator compositePropensityCreator;

    public CooperativeZPropensityCreator(CompositePropensityCreator compositePropensityCreator) {
        super(InfoType.COOP_Z_PROPENSITY);
        this.compositePropensityCreator = compositePropensityCreator;
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        if (term != null) return term;
        /* load parameters if this is the first time the function is called */
        if (parameters == null)
            parameters = new CooperativeZPropensityParameters(commands);
        CompositePropensityEnergy compositePropensityEnergy = (CompositePropensityEnergy) compositePropensityCreator.term();
        if (compositePropensityEnergy == null)
            throw new RuntimeException("Apparently an attempt to insantiate CooperativeZPropensity " +
                    "before CompositePropensityEnergy");
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.COOP_Z_PROPENSITY, " A cooperative term that effect the distribution of propensity values.", weight);
        return term = new CooperativeZPropensityEnergy(compositePropensityEnergy, energyInfoElement, parameters);
    }

}


