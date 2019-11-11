package meshi.energy.contacts;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.util.filters.HydrophobicSideChains;
import meshi.util.filters.PolarSideChains;
import meshi.util.filters.SecondaryStructureFilter;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.CommandList;
import meshi.util.filters.Filter;
import meshi.util.info.InfoType;

/**
 * Created by chen on 28/02/2015.
 */
public class ContactsCreator extends EnergyCreator {
    public ContactsCreator()  {
        super(InfoType.CONTACTS_ETC);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        Filter hydrophobicSideChains = new HydrophobicSideChains(false);
        Filter secondaryStructureFilter = new SecondaryStructureFilter();
        Filter polarFilter = new PolarSideChains();

        AtomList nonPolarSS = new AtomList(protein.molecularSystem);
        AtomList polarSS    = new AtomList(protein.molecularSystem);
        for (int iChain = 0; iChain < protein.chains().size(); iChain++) {
            nonPolarSS.addAll(protein.chains().get(iChain).atoms().filter(hydrophobicSideChains).filter(secondaryStructureFilter));
            polarSS.addAll(protein.chains().get(iChain).atoms().filter(polarFilter).filter(secondaryStructureFilter));
        }
        term = new ContactsAndSASA(polarSS,nonPolarSS,new ContactsInfo(infoType,
                                                                       "Course grain evaluation of hydrophobic hydrophilic relations"));
        return term;
    }
}
