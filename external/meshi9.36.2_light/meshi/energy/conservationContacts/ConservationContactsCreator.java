package meshi.energy.conservationContacts;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 */
public class ConservationContactsCreator extends EnergyCreator {
    private boolean fileFound = false;
    public ConservationContactsCreator(InfoType infoType) {
        super(infoType);
    }

    public  void setFileFound(boolean fileFound) {
        this.fileFound = fileFound;
    }
      public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands){
          return new ConservationContactsRatio(protein.atoms().CAFilter(),
                                               new EnergyInfoElement(infoType,
                                               "Ratio between the number of contacts by conserved and other residues.",weight),fileFound);
      }
}
