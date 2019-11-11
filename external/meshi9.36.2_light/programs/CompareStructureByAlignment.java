package programs;

import meshi.PDB.PdbLineATOM;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.ca.CaResidue;
import meshi.molecularElements.extendedAtoms.BackboneResidue;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.sequences.SequenceList;
import meshi.util.Rms;
import meshi.util.filters.Filter;

import java.io.File;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 06/06/2010
 * Time: 14:44:29
 * .
 */
public class CompareStructureByAlignment {
    public static final String NAME = "CompareStructureByAlignment";

    public static void main(String[] args) {
        String alignmentFile = args[0];      //fasta file with the 2 sequences that are aligned after deleting the gaps
        String dirName1 = args[1];
        File dir1 = new File(dirName1);      // directory of the models of the first target
        String dirName2 = args[2];
        File dir2 = new File(dirName2);    // directory of the models of the second target
        SequenceList sequenceList = new SequenceList(alignmentFile);
        int sequenceAligmentLength = 120; //Todo change it to a parameter - this is the number of the residues that are part of the alignment - any protein that is shorter the this number should be removed.


        if (!dir1.exists() || !dir2.exists()) throw new RuntimeException(dir1 + " or " + dir2 + " does not exist.");
        if (!dir1.isDirectory() || !dir2.isDirectory())
            throw new RuntimeException(dir1 + " or " + dir2 + " is not a directory.");
        File[] filesDir1 = dir1.listFiles();
        File[] filesDir2 = dir2.listFiles();
        boolean[] excluded1 = new boolean[filesDir1.length];
        boolean[] excluded2 = new boolean[filesDir2.length];
        Protein model1;
        Protein model2;
        ResidueAlignment residueAlignment;

        double[][] rmsValues = new double[filesDir1.length][filesDir2.length];
        double[][] gdtValues = new double[filesDir1.length][filesDir2.length];
        ArrayList<String> ans = new ArrayList<String>();

        double bestRms = 10000;
        double bestGdt = -1;
        Protein bestRmsModel1 = null, bestRmsModel2 = null, bestGdtModel1 = null, bestGdtModel2 = null;


        test(filesDir1, sequenceAligmentLength);
        test(filesDir2, sequenceAligmentLength);
        for (int i = 0; i < filesDir1.length; i++) {
            if (bestRmsModel1 != null)
                System.out.println("Best Rms = " + bestRms + " " + bestRmsModel1.name() + " " + bestRmsModel2.name());
            if (bestGdtModel1 != null)
                System.out.println("Best GDT_TS = " + bestGdt + " " + bestGdtModel1.name() + " " + bestGdtModel2.name());
            if (filesDir1[i] != null) {
                System.out.println("analyzing  " + filesDir1[i].getName() + " " + i + " of " + filesDir1.length);
                try {

                    new MolecularSystem();
                    model1 = new Protein(filesDir1[i].getAbsolutePath(), new PdbLineATOM(), CaResidue.creator, null);
                    for (int j = 0; j < filesDir2.length; j++) {
                        if (filesDir2[j] != null) {
                            System.out.print(".");
                            if (j % 50 == 0)
                                System.out.print(j);
                            if (j % 300 == 0) System.out.println();
                            try {
                                new MolecularSystem();
                                model2 = new Protein(filesDir2[j].getAbsolutePath(), new PdbLineATOM(), CaResidue.creator, null);

                                //********************************************************************************************************************************************************************************************
                                /************************************************************* this is the section that deal with the rms calculation ****************************************
                                 /**************************************************************************************************************************************************************************************************/
                                residueAlignment = new ResidueAlignment(model1.chain(), model2.chain(), sequenceList);
                                double gdtValue = -1;
                                Rms rms = new Rms(residueAlignment);
                                double rmsValue = rms.getRms();
                                try {
                                    gdtValue = Rms.gdt(residueAlignment, new IsCa(), residueAlignment.size())[0];
                                }
                                catch (Exception ex) {
                                    for (Atom atom : model2.chain().atoms()) {
                                        System.out.println("xxxx " + atom + "; " + atom.molecularSystem);
                                    }
                                    residueAlignment.print();
                                    throw new RuntimeException(ex);
                                }
                                rmsValues[i][j] = rmsValue;
                                gdtValues[i][j] = gdtValue;
                                if (rmsValue < bestRms) {
                                    bestRms = rmsValue;
                                    bestRmsModel1 = model1;
                                    bestRmsModel2 = model2;
                                }
                                if (gdtValue > bestGdt) {
                                    bestGdt = gdtValue;
                                    bestGdtModel1 = model1;
                                    bestGdtModel2 = model2;
                                }
                                if (rmsValue < 8 && gdtValue > 0.5) {
                                    ans.add(model1.name() + " " + model2.name() + " rms: " + rmsValue + " gdt: " + gdtValue);
                                }
                            }
                            catch (Exception ex) {
                                System.out.println("a probelm  with " + filesDir1[i].getName() + " " + filesDir2[j].getName() + " thus it was not considered!");
                                throw new RuntimeException(ex);

                            }
                        }
                        model2 = null;
                    }
                    System.out.println();

                }
                catch (Exception ex) {
                    System.out.println("a probelm  with " + filesDir1[i].getName() + " thus it was not considered!");
                    ex.printStackTrace();
                    throw new RuntimeException(ex);
                }
                model1 = null;
            }
        }

        for (String s : ans) {
            System.out.println(s);
        }
    }

    public static void test(File[] filesDir1, int minimalLength) {
        for (int i = 0; i < filesDir1.length; i++) {
            System.out.println("Testing  " + filesDir1[i].getName() + " " + i + " of " + filesDir1.length);
            Protein model1;
            try {
                new MolecularSystem();
                model1 = new Protein(filesDir1[i].getAbsolutePath(), new PdbLineATOM(), BackboneResidue.creator, null);
                ResidueAlignment residueAlignment = new ResidueAlignment(model1.chain(), "xxx", model1.chain(), "xxx", ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                double[] gdtValue = Rms.gdt(residueAlignment, new IsCa(), residueAlignment.size());
            }
            catch (Exception ex) {
                System.out.println("Failed to open create a protein from " + filesDir1[i].getAbsolutePath() + " due to " + ex);
                filesDir1[i] = null;
                continue;
            }
            if (model1.chain().numberOfNonDummyResidues() < minimalLength) {
                System.out.println("dir1: protein " + model1.name() + " is too short: " + model1.chain().numberOfNonDummyResidues());
                filesDir1[i] = null;
            }
        }
    }


    private static class IsBackbone implements Filter {
        public boolean accept(Object obj) {
            if (obj == null) return false;
            Atom tmp = (Atom) obj;
            return tmp.isBackbone();
        }
    }

    private static class IsCa implements Filter {
        public boolean accept(Object obj) {
            if (obj == null) return false;
            Atom tmp = (Atom) obj;
            return tmp.backboneCA();
        }
    }
}
