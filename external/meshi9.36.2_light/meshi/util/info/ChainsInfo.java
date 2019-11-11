package meshi.util.info;

import meshi.energy.AbstractEnergy;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.util.ResidueData;

import java.util.ArrayList;

/**
 * Created by chen on 20/08/2017.
 */
public class ChainsInfo extends ArrayList<ChainInfo> {
    public ChainsInfo() {}

    public ChainsInfo(Protein protein) {
        for (Chain chain : protein.chains()) {
            add(new ChainInfo(chain));
        }
    }

    public void add(ResidueData residueData, InfoType type) {
        for (ChainInfo chainInfo : this) {
            for (ResidueInfo residueInfo : chainInfo) {
                if (! residueInfo.dummy()) {
                    residueInfo.add(new DoubleInfoElement(type, type.tag, residueData.get(residueInfo.chainNumber(), residueInfo.number())));
                }
            }
        }
    }
}
