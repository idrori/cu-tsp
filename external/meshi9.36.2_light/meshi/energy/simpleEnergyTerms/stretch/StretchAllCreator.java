package meshi.energy.simpleEnergyTerms.stretch;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 30/04/14
 * Time: 23:21
 * To change this template use File | Settings | File Templates.
 */
public class StretchAllCreator extends EnergyCreator{
    public StretchAllCreator() {
        super(InfoType.STRETCH);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix dm, CommandList commands){

        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.STRETCH_ALL, "stretchAll", weight);
        term =  new StretchAll(protein,energyInfoElement);
        return term;

    }
}
