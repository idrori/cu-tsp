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
 * Date: 06/08/2009
 * Time: 12:00:30
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZStdPropensityCreator extends EnergyCreator implements KeyWords {

    private static CooperativeZPropensityParameters parameters = null;
    private CompositePropensityCreator propensityEnergyCreator;
    private CooperativeZPropensityCreator cooperativeZPropensityCreator;

    public CooperativeZStdPropensityCreator(CooperativeZPropensityCreator cooperativeZPropensityCreator) {
        super(InfoType.COOP_Z_STD_PROPENSITY);
        this.cooperativeZPropensityCreator = cooperativeZPropensityCreator;
        this.propensityEnergyCreator = cooperativeZPropensityCreator.compositePropensityCreator;
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        if (term != null) return term;
        /* load parameters if this is the first time the function is called */
        parameters = cooperativeZPropensityCreator.parameters;
        if (parameters == null)
            throw new RuntimeException("No parameters for " + this);


        CompositePropensityEnergy propencityEnergy = (CompositePropensityEnergy) propensityEnergyCreator.term();
        if (propencityEnergy == null)
            throw new RuntimeException("Apparently an attempt to insantiate CooperativeZstdPropensity " +
                    "before CompositePropensityEnergy");
        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.COOP_Z_STD_PROPENSITY, " A cooperative term that effect the distribution of propensity values.", weight);
        return term = new CooperativeZStdPropensityEnergy(propencityEnergy, energyInfoElement, parameters);
    }

}
