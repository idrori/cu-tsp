package meshi.energy.one;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.atoms.Atom;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 13/12/11
 * Time: 16:40
 * To change this template use File | Settings | File Templates.
 */
public class One extends AbstractEnergy{

    public One(EnergyInfoElement info) {
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
    }
    public EnergyInfoElement evaluate() {
        info.setValue(1.0);
        return  info;
    }
     public void evaluateAtoms() {}
    public void test(TotalEnergy totalEnergy, Atom atom){
        System.out.println("Testing One");
    }
}
