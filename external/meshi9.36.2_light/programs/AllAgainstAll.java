package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.MeshiProgram;
import meshi.util.Rms;
import meshi.util.UpdateableException;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;

/**
 * Created by chen on 21/04/2015.
 */
public class AllAgainstAll {
    public static void main(String[] args) throws Exception{
        File currentDir = new File(".");
        File[] files = currentDir.listFiles();
        for (File file : files)
            if (file.isDirectory()) {
                MeshiWriter writer = new MeshiWriter(file.getName()+"_similarityTable.dat");
                analyzeTarget(file,writer);
                writer.close();
            }
    }

    public static void analyzeTarget(File targetDir, MeshiWriter writer) throws Exception{
        Protein iProtein, jProtein;
        File[] files = targetDir.listFiles();
        for (int iFile = 0; iFile < files.length; iFile++) {
            String iFileName = files[iFile].getAbsolutePath();
            if (iFileName.endsWith(".pdb")) {
                iProtein = new Protein(new AtomList(iFileName), new ResidueExtendedAtomsCreator());
                for (int jFile = 0; jFile < files.length; jFile++) {
                    String jFileName = files[jFile].getAbsolutePath();
                    if (jFileName.endsWith(".pdb")) {
                        try {
                            jProtein = new Protein(new AtomList(jFileName), new ResidueExtendedAtomsCreator());
                        } catch (Exception ex) {
                            System.out.println("Problem while building a model from " + files[jFile]);
                            ex.printStackTrace();
                            throw (ex);
                        }
                        double rms = Rms.rms(iProtein, jProtein, ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                        double[] gdt = Rms.gdt(iProtein, jProtein);
                        writer.println(iProtein.name() + "\t" + jProtein.name() + "\t" + rms + "\t" + gdt[0]);
                        System.out.print(".");
                    }
                }
                System.out.println();
                System.out.print(iFileName+"  ");
            }
        }
    }
}
