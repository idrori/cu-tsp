package meshi.util.info;

import meshi.energy.TotalEnergy;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;

import java.util.ArrayList;

/**
 * Created by chen on 17/08/2017.
 */
public class ChainInfo extends ArrayList<ResidueInfo> {
    public ChainInfo(Chain chain) {
        super();
        for (Residue residue : chain) {
            add(new ResidueInfo(residue));
        }
    }
}
