package meshi.energy.contacts;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

import java.util.ArrayList;

/**
 * Created by chen on 29/06/2015.
 */
public class ContactsInfo extends EnergyInfoElement {
    public final MeshiInfo contacts12, contacts14, contacts14CoreRatio, sasaRatio;

    public ContactsInfo(InfoType infoType, String comment) {
        super(infoType, comment, 1);
        ArrayList<MeshiInfo> meshiInfos = getChildren();
        meshiInfos.add(contacts12 = new MeshiInfo(InfoType.CONTACTS12,
                "Z-score of length normalized average number of 12A contacts between hydrophobic side-chains in SS"));
        meshiInfos.add(contacts14 = new MeshiInfo(InfoType.CONTACTS14,
                "Z-score of length normalized average number of 14A contacts between hydrophobic side-chains in SS"));
        meshiInfos.add(contacts14CoreRatio = new MeshiInfo(InfoType.CONTACTS14_CORE_RATIO,
                "Mahalanobis distance of length normalized average number of 14A contacts of all hydrophobic side-chains in SS "+
                        " and the average number of contacts among the most compact getAtoms (core)."));
        meshiInfos.add(sasaRatio = new MeshiInfo(InfoType.SASA_RATIO,
                "Mahalanobis distance of length normalized average SASA (as measured by DSSP)of all hydrophobic side-chains in SS "+
                        " and the polar ones."));
    }
}
