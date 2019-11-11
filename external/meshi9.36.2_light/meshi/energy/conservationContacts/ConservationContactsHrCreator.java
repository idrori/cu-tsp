package meshi.energy.conservationContacts;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.CommandList;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

/**
 */
public class ConservationContactsHrCreator extends EnergyCreator {
    boolean fileFound = false;
    public ConservationContactsHrCreator() {
        super(InfoType.CONSERVATION_CONTACTS_HR);
    }

    public  void setFileFound(boolean fileFound) {
        this.fileFound = fileFound;
    }
      public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands){
          EnergyInfoElement contactsInfo;
          EnergyInfoElement rgRatioInfo;
          MeshiInfo[] meshiInfos = new MeshiInfo[2];
//          contactsInfo = new EnergyInfoElement(InfoType.CONTACTS_HR,"Average number of CA contacts (8A from a CA of a residue that is more than 10 residues away");
//          rgRatioInfo  = new EnergyInfoElement(InfoType.CONSERVATION_RG_RATIO_HR,"The ratio between the RGs of conserved CAs and all CAs");
//          meshiInfos[0] = contactsInfo;
//          meshiInfos[1] = rgRatioInfo;

          return new ConservationContactsHrRatio(protein.atoms().sideChains(),
                                                 new EnergyInfoElement(infoType, "Ratio between the number of contacts by conserved and other residues.",weight),
                                                 fileFound);
      }
}
