package meshi.applications.conserv.statistic;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

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
          return new ConservationContactsHrRatio(protein.atoms().sideChains(),
        		  new EnergyInfoElement(infoType,
                        "Ratio between the number of contacts by conserved and other residues.",weight),fileFound);
      }
}
