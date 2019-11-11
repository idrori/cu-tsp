package meshi.energy.conservationContacts;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Conservation;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.Utils;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;

import java.util.ArrayList;

/**

 */
public class ConservationContactsRatio extends AbstractEnergy {
    AtomList atoms,hAtoms;
    EnergyInfoElement contactsInfo;
    EnergyInfoElement rgRatioInfo;
    EnergyInfoElement hRgRatioInfo;
    private double contactDistance;
    private InfoType conservationRgRatio,conservationRgRatioH,contactsType;
    int[] contacts;
    private boolean fileFound;
    ArrayList<Double> weights,hWeights;
    int nAtoms;
    public ConservationContactsRatio(AtomList atomList, EnergyInfoElement info, boolean fileFound) {
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
        conservationRgRatio = InfoType.CONSERVATION_RG_RATIO;
        conservationRgRatioH = InfoType.CONSERVATION_H_RG_RATIO;
        switch (info.type){
            case CONSERVATION_CONTACTS8:
                contactDistance = 8;
                contactsType = InfoType.CONTACTS8;
                break;
            case CONSERVATION_CONTACTS11:
                contactDistance = 11;
                contactsType = InfoType.CONTACTS11;
                break;
            case CONSERVATION_CONTACTS15:
                contactsType = InfoType.CONTACTS15;
                contactDistance = 15;
                break;
            default: throw new RuntimeException("This is weird."+info.type);
        }
        atoms = atomList;
        contactsInfo = new EnergyInfoElement(contactsType,"Average number of CA contacts ("+contactDistance+" from a CA of a residue that is more than 10 residues away");
        rgRatioInfo = new EnergyInfoElement(conservationRgRatio,"The ratio between the RGs of conserved CAs and all CAs");
        hRgRatioInfo = new EnergyInfoElement(conservationRgRatioH,"The ratio between the RGs of conserved CAs and all CAs - for hydrophobic residues");
        info.getChildren().add(contactsInfo);
        if (contactDistance == 8) {
            info.getChildren().add(rgRatioInfo);
            info.getChildren().add(hRgRatioInfo);
        }
        this.fileFound = fileFound;
        comment = "conserfEnergy";
        nAtoms = atoms.size();
        contacts = new int[nAtoms];
        hAtoms = new AtomList(atoms.molecularSystem);
        weights = new ArrayList<Double>();
        hWeights = new ArrayList<Double>();
        if (fileFound) {
           for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                Atom atomI = atoms.atomAt(iAtom);
                Conservation conservation = atomI.residue().conservation();
                if (conservation == null) throw new RuntimeException("Conservation is null in: "+atomI);
                double weight = conservation.conservationWeight();
                weights.add(new Double(weight));
               if (atomI.residue().type.isHydrophobic) {
                   hAtoms.add(atomI);
                   hWeights.add(new Double(weight));
               }
           }
        }
    }

    public EnergyInfoElement evaluate() {
        double  countAll = 0, countWeighted = 0, sumWeights= 0;
        for (int iAtom = 0; iAtom < nAtoms; iAtom++)
            contacts[iAtom] = 0;
        for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
            Atom atomI = atoms.atomAt(iAtom);
            if (atomI.nowhere()) continue;
            for (int jAtom = iAtom+10; jAtom< nAtoms; jAtom++) {
                Atom atomJ = atoms.atomAt(jAtom);
                if (atomJ.nowhere()) continue;
                double d = (new FreeDistance(atomI,atomJ)).distance();
                 if (d < contactDistance) {
                     contacts[iAtom]++;
                     contacts[jAtom]++;
                     countAll += 2;
                 }
            }
        }

        double rg, hRg,rgRatio,hRgRatio,energy;
        double weightedRg,hWeightedRg;
        if (fileFound) {
           for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                double weight = weights.get(iAtom).doubleValue();
                countWeighted += contacts[iAtom]*weight;
                sumWeights += weight;
           }

           rg = Utils.radiusOfGyration(atoms);
           weightedRg = Utils.radiusOfGyration(atoms,weights);
           hRg = Utils.radiusOfGyration(hAtoms);
           hWeightedRg = Utils.radiusOfGyration(hAtoms,hWeights);
           rgRatio = 1.0*weightedRg/rg;
           hRgRatio = 1.0*hWeightedRg/hRg;
           energy = (1.0*countAll/nAtoms)/(countWeighted/sumWeights);
           info.setValue(new Double(energy));
        }
        else {
            info.setValue(0.0);
            rgRatioInfo.setValue(0.0);
            hRgRatioInfo.setValue(0.0);
        }
        contactsInfo.setValue(new Double(1.0*countAll / nAtoms));
        return info;
    }


    public void evaluateAtoms() {}
    public void test(TotalEnergy energy, Atom atom) {}
}
