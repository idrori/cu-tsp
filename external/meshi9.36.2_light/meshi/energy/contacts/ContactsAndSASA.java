package meshi.energy.contacts;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.CommandList;
import meshi.util.Utils;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by chen on 28/02/2015.
 */
public class ContactsAndSASA extends AbstractEnergy {
    private static final double[][] CONTACTS12_LINEAR = {{-9.608,-6.006}};
    private static final double[] CONTACTS12_MU = {0.99944};
    private static final double[][] CONTACTS12_INV_SIGMA = {{54.4194}};

    private static final double[][] CONTACTS14_LINEAR = {{-16.9879,14.2102}};
    private static final double[] CONTACTS14_MU = {0.99953};
    private static final double[][] CONTACTS14_INV_SIGMA = {{91.6888}};


    private static final double[][] CONTACTS14_CORE_TOTAL_LINEAR = {{-23.4241, 27.8281},{-16.3667, 10.7837}};
    private static final double[] CONTACTS14_CORE_TOTAL_MU = {0.99951, 0.9995};
    private static final double[][] CONTACTS14_CORE_TOTAL_INV_SIGMA = {{2764.213, -2392.8542},{-2392.8542, 2164.6184}};


    private static final double[][] SASA_LINEAR = {{-0.052308, 0.44021},{-0.072959, 0.49148}};
    private static final double[] SASA_MU = {0.99981, 0.99991};
    private static final double[][] SASA_INV_SIGMA = {{76.1307, -19.7917},{-19.7917, 32.7084}};



    private final AtomList polarAtoms, nonPolarAtoms;
    private Mahalanobis contacts12 = new Mahalanobis(CONTACTS12_LINEAR,CONTACTS12_MU,CONTACTS12_INV_SIGMA);
    private Mahalanobis contacts14 = new Mahalanobis(CONTACTS14_LINEAR, CONTACTS14_MU, CONTACTS14_INV_SIGMA);
    private Mahalanobis contacts14CoreTotal = new Mahalanobis(CONTACTS14_CORE_TOTAL_LINEAR,CONTACTS14_CORE_TOTAL_MU,CONTACTS14_CORE_TOTAL_INV_SIGMA);
    private Mahalanobis SASAratio = new Mahalanobis(SASA_LINEAR,SASA_MU, SASA_INV_SIGMA);


    public ContactsAndSASA(AtomList polarAtoms, AtomList nonPolarAtoms, EnergyInfoElement info) {
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
        this.polarAtoms = polarAtoms;
        this.nonPolarAtoms = nonPolarAtoms;
        comment = "contactsAndSASA";
        ////////////////contacts12 = new Mahalanobis()
    }

    //For debugging
    public String normalizers() {
       return ""+Math.log(nonPolarAtoms.size()+1)+" "+Math.log(polarAtoms.size()+1);
    }
    public static double[] getContacts(AtomList atoms, double threshold) {
        int[] nContacts = new int[atoms.size()];
        for (int i = 0; i < atoms.size(); i++)  {
            Atom atomI = atoms.get(i);
            for (int j = i+1; j < atoms.size();j++) {
                Atom atomJ = atoms.get(j);
                if (atomI.distanceFrom(atomJ) <= threshold) {
                    nContacts[i]--;
                    nContacts[j]--;
                }
            }
        }
        double[] out = new double[2];
        double sum = 0;
        Arrays.sort(nContacts);

        for (int i = 0; i < atoms.size()/2; i++)
            sum += nContacts[i];
        out[0] = sum/((atoms.size()+1.0)/2);
        for (int i = atoms.size()/2; i < atoms.size(); i++)
            sum += nContacts[i];
        out[1] = sum/(atoms.size()+1.0);

        return out;
    }

    public static double[] getAccessibility(AtomList atoms) {
        double[] out = new double[2];
        Residue residue = null;
        ArrayList<Double> accessibility = new ArrayList();
        int iResidue = 0;
        for (int i = 0; i < atoms.size(); i++)  {
            Atom atomI = atoms.get(i);
            if (atomI.residue() != residue) {
                residue = atomI.residue();
                iResidue++;
                boolean found = false;
                for (Atom atom : residue.getAtoms())
                    if ((!atom.isHydrogen()) && atom.nowhere()) found = true;
                if (found) accessibility.add(1.0); // A missing atom indicates a highly exposed residue but dssp do not count the getAccessibility of missing atoms.
                else accessibility.add(residue.getRelativeAccesibility());
            }
         }
        double sum = 0;
        Double[] accessibilityArray = accessibility.toArray(new Double[accessibility.size()]);
        Arrays.sort(accessibilityArray);

        for (int i = 0; i < iResidue/2; i++)
            sum += accessibilityArray[i];
        out[0] = sum/(iResidue/2+1);
        for (int i = iResidue/2; i < iResidue; i++)
            sum += accessibilityArray[i];
        out[1] = sum/(iResidue+1);

        return out;
    }

    public void evaluateAtoms() {}

    public EnergyInfoElement evaluate() {
        double logNnonPolars = Math.log(nonPolarAtoms.size()+1);
        double logNpolar     = Math.log(polarAtoms.size()+1);
        if (logNnonPolars == 0)
            logNnonPolars = 0.001;
        if (logNpolar == 0)
            logNpolar = 0.001;

        double[] normalizers1 = new double[1];
        double[] normalizers2 = new double[2];
        double[] features1 = new double[1];
        double[] features2 = new double[2];

        double[] hydrophobicSideChainsSsMContacts12 = getContacts(nonPolarAtoms,12);
        double[] hydrophobicSideChainsSsMContacts14 = ContactsAndSASA.getContacts(nonPolarAtoms, 14);
        double[] accesiblityHhydrophobicSs = getAccessibility(nonPolarAtoms);
        double[] accesiblityPolarSs = getAccessibility(polarAtoms);

        normalizers1[0] = logNnonPolars;
        features1[0]    = hydrophobicSideChainsSsMContacts12[1];

        double contacts12Value = contacts12.distance(normalizers1,features1);

        features1[0]    = hydrophobicSideChainsSsMContacts14[1];
        double contacts14Value = contacts14.distance(normalizers1,features1);

        normalizers2[0] = logNnonPolars;
        normalizers2[1] = logNnonPolars;
        features2[0]    = hydrophobicSideChainsSsMContacts14[0];
        features2[1]    = hydrophobicSideChainsSsMContacts14[1];
        double contacts14coreTotalValue = contacts14CoreTotal.distance(normalizers2,features2);

        normalizers2[1] = logNpolar;
        features2[0] = accesiblityHhydrophobicSs[1];
        features2[1] = accesiblityPolarSs[0];
        double SASA = SASAratio.distance(normalizers2,features2);



        double energy = -contacts12Value-0.5*contacts14Value-0.5*contacts14coreTotalValue+0.5*SASA;
        if (Double.isNaN(energy)) energy = -9999.999;
        info.setValue(energy);
        ((ContactsInfo) info).contacts12.setValue(contacts12Value);
        ((ContactsInfo) info).contacts14.setValue(contacts14Value);
        ((ContactsInfo) info).contacts14CoreRatio.setValue(contacts14coreTotalValue);
        ((ContactsInfo) info).sasaRatio.setValue(SASA);
        return info;

    }
    public void test(TotalEnergy totalEnergy, Atom atom) {}




    private static class Mahalanobis {
        private double [][]linear;
        private double [] mu;
        private double [][] invSigma;

        private Mahalanobis(double[][] linear, double[] mu, double[][] invSigma) {
            this.linear = linear;
            this.mu = mu;
            this.invSigma = invSigma;
            if(mu.length  != invSigma.length)   throw new RuntimeException("This is weird");
            if(mu.length  != invSigma[0].length)   throw new RuntimeException("This is weird");
            if (linear != null) {
                if (linear.length > mu.length) throw new RuntimeException("This is weird");
                if (linear[0].length != 2) throw new RuntimeException("This is weird");
            }
        }


        public double distance(double normalizers[], double[] features) {
            int nNormalizers;
            int nFeatures = features.length;
            double[] transformed = new double[features.length];
            double[] temp = new double[features.length];
            for (int i = 0; i < features.length; i++)
                temp[i] = 0;
            if (nFeatures != mu.length) throw new RuntimeException("This is weird");
            if (normalizers != null) {
                nNormalizers = normalizers.length;
                if (nNormalizers > nFeatures) throw new RuntimeException("This is weird");
                for (int iFeature = 0; iFeature < normalizers.length; iFeature++) {
                        transformed[iFeature] = (features[iFeature] - linear[iFeature][1]) / normalizers[iFeature] / linear[iFeature][0];
               }
            } else {
                nNormalizers = 0;
            }

            for (int iFeature = nNormalizers; iFeature < features.length; iFeature++) {
                    transformed[iFeature] = features[iFeature];
            }

            for (int iFeature = 0; iFeature < features.length; iFeature++) {
                    transformed[iFeature] = transformed[iFeature] - mu[iFeature];
            }

            double distance = 0;

                for (int iFeature = 0; iFeature < features.length; iFeature++) {
                    for (int jFeature = 0; jFeature < features.length; jFeature++) {
                        temp[iFeature] += transformed[jFeature] * invSigma[iFeature][jFeature];
                    }
                }
            for (int iFeature = 0; iFeature < features.length; iFeature++)
                distance += temp[iFeature]*transformed[iFeature];

            return Math.sqrt(distance);
        }

    }

    public static void main(String[] args) {
        new MolecularSystem();
        Protein protein = new Protein(new AtomList("T0388-D1.3Dpro_TS4.pdb"), new ResidueExtendedAtomsCreator());
        Utils.AssignDSSP(protein,"T0388-D1.3Dpro_TS4.pdb.dssp" );
        ContactsCreator creator = new ContactsCreator();
        ContactsAndSASA energyTerm = (ContactsAndSASA) creator.createEnergyTerm(protein,null,new CommandList("commands.MCM"));
        System.out.println(energyTerm.normalizers());
        EnergyInfoElement info = energyTerm.evaluate();
        System.out.println(info.toXml());
    }
}
