package meshi.energy.one;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 13/12/11
 * Time: 21:13
 * To change this template use File | Settings | File Templates.
 */
public class OneCreator extends EnergyCreator {
    private EnergyInfoElement info = new EnergyInfoElement(InfoType.ONE,"The constant 1. Useful for regression",weight);
    public OneCreator() {
        super(InfoType.ONE);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix dm, CommandList cl) {
        return new One(info);
    }

}
