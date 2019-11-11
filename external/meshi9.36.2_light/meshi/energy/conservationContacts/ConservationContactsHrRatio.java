package meshi.energy.conservationContacts;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

/**

 */
public class ConservationContactsHrRatio extends AbstractEnergy {
    AtomList atoms;
    EnergyInfoElement contactsInfo;
    EnergyInfoElement rgRatioInfo;
    MeshiInfo[] meshiInfos = new MeshiInfo[2];
    private boolean fileFound;
    public ConservationContactsHrRatio(AtomList atomList, EnergyInfoElement info, boolean fileFound) {
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
        atoms = atomList;
        contactsInfo = new EnergyInfoElement(InfoType.CONTACTS_HR,"Average number of CA contacts (8A from a CA of a residue that is more than 10 residues away");
        rgRatioInfo  = new EnergyInfoElement(InfoType.CONSERVATION_RG_RATIO_HR,"The ratio between the RGs of conserved CAs and all CAs");
        info.getChildren().add(contactsInfo);
        info.getChildren().add(rgRatioInfo);
        this.fileFound = fileFound;
        comment = "conserfHrEnergy";
    }

    public EnergyInfoElement evaluate() {
        double  countAll = 0, countWeighted = 0, sumWeights= 0;
        int N = atoms.size();
        int[] contacts = new int[N];
        for (int iAtom = 0; iAtom < N; iAtom++) {
            Atom atomI = atoms.atomAt(iAtom);
            if (!atomI.nowhere()) {
                Residue residueI = atomI.residue();
                for (int jAtom = iAtom + 1; jAtom < N; jAtom++) {
                    Atom atomJ = atoms.atomAt(jAtom);
                    if (!atomJ.nowhere()) {
                        if (atomJ.residue() != residueI) {
                            double d = (new FreeDistance(atomI, atomJ)).distance();
                            if (d < 6) {
                                contacts[iAtom]++;
                                contacts[jAtom]++;
                            }
                        }
                    }
                }
            }
        }
        double xCmWeighted = 0, yCmWeighted = 0, zCmWeighted = 0;
        double xCm = 0, yCm = 0, zCm = 0;
        for (int iAtom = 0; iAtom < N; iAtom++) {
                Atom atomI = atoms.atomAt(iAtom);
                countAll += contacts[iAtom];
                if (fileFound) {
                    double weight = atomI.residue().conservation().conservationWeight();
                    countWeighted += contacts[iAtom]*weight;
                    sumWeights += weight;
                    xCmWeighted += atomI.x()*weight;
                    yCmWeighted += atomI.y()*weight;
                    zCmWeighted += atomI.z()*weight;
                    xCm += atomI.x();
                    yCm += atomI.y();
                    zCm += atomI.z();
                }
         }
        if(fileFound) {
            xCmWeighted /= sumWeights;
            yCmWeighted /= sumWeights;
            zCmWeighted /= sumWeights;
            xCm /= N;
            yCm /= N;
            zCm /= N;
            double sum = 0, weightedSum = 0;
            for (int iAtom = 0; iAtom < N; iAtom++) {
                Atom atomI = atoms.atomAt(iAtom);
                double weight = atomI.residue().conservation().conservationWeight();
                double x = atomI.x();
                double y = atomI.y();
                double z = atomI.z();
                sum +=  (xCm-x)*(xCm-x);
                sum +=  (yCm-y)*(yCm-y);
                sum +=  (zCm-z)*(zCm-z);
                weightedSum +=  (xCmWeighted -x)*(xCmWeighted -x)*weight;
                weightedSum +=  (yCmWeighted -y)*(yCmWeighted -y)*weight;
                weightedSum +=  (zCmWeighted -z)*(zCmWeighted -z)*weight;
            }
            double rg = Math.sqrt(sum/N);
            double weightedRg = Math.sqrt(weightedSum/sumWeights);
            double rgRatio = weightedRg/rg;
            double energy = (1.0*countAll/N)/(countWeighted/sumWeights);
            info.setValue(energy);
            rgRatioInfo.setValue(rgRatio);
        }
        else {
            info.setValue(0.0);
            rgRatioInfo.setValue(0.0);
        }
        contactsInfo.setValue(1.0*countAll / N);
        return info;
    }


    public void evaluateAtoms() {}
    public void test(TotalEnergy energy, Atom atom) {}
}
