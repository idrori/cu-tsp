package meshi.applications.conserv.statistic;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.parameters.AtomType;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

/**

 */
public class ConservationContactsHrRatio extends AbstractEnergy {

	AtomList 			atoms;
	MeshiInfo contactsInfo;
	MeshiInfo 	rgRatioInfo;
	private boolean 	fileFound;
	private double 		mRadius;

	public ConservationContactsHrRatio(AtomList atomList, EnergyInfoElement info, boolean fileFound) {
		super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
		mRadius = 6;
		atoms = atomList;
		contactsInfo = new MeshiInfo(InfoType.CONTACTS_HR, Double.valueOf(-9999.999),"Average number of CA contacts (8A from a CA of a residue that is more than 10 residues away");
		rgRatioInfo  = new MeshiInfo(InfoType.CONSERVATION_RG_RATIO_HR,  Double.valueOf(-9999.999),"The ratio between the RGs of conserved CAs and all CAs");
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
				for (int jAtom = iAtom+1; jAtom< N; jAtom++) {
					Atom atomJ = atoms.atomAt(jAtom);
					if ((atomJ.residue() != residueI) && (!atomJ.nowhere()))  {
						double d = (new FreeDistance(atomI,atomJ)).distance();
						if (d < mRadius) {
							contacts[iAtom]++;
							contacts[jAtom]++;
						}
					}
				}
			}
			else if ((atomI.type() != AtomType.TRO) ){
                throw new RuntimeException(
                        "!!!!!! Chen do something please !!!!!!!\n"+
                        atomI+"\n"+
                        atomI.type());
            }
		}
		double xCmWeighted = 0, yCmWeighted = 0, zCmWeighted = 0;
		double xCm = 0, yCm = 0, zCm = 0;
		for (int iAtom = 0; iAtom < N; iAtom++) {
			Atom atomI = atoms.atomAt(iAtom);
			countAll += contacts[iAtom];
			if (fileFound & (!atomI.nowhere())) {
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
				if (!atomI.nowhere()) {
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
			}
			double rg = Math.sqrt(sum/N);
			double weightedRg = Math.sqrt(weightedSum/sumWeights);
			double rgRatio = weightedRg/rg;
			double energy = (countAll/N)/(countWeighted/sumWeights);
			if (energy > 10000) {
				System.out.println("AAAAA");
			}
			info.setValue(energy);
			rgRatioInfo.setValue(rgRatio);
		}
		else {
			info.setValue(0);
			rgRatioInfo.setValue(Double.MAX_VALUE);
		}
		contactsInfo.setValue(countAll / N);
		System.out.println("Number of getAtoms is : " + N + "   Contacts: " + contactsInfo.getValue() + "   Rg: " + rgRatioInfo.getValue());
		return info;
	}

	/**
	 * Set the radius for calculating, if this method is never called, defaulrt radius is 6.
	 * @param radius The max radius for calculating
	 */
	public void setRadius(double radius) {
		mRadius = radius;
	}


	public void evaluateAtoms() {}
	public void test(TotalEnergy energy, Atom atom) {}
}
