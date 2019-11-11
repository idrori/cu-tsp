package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.Rms;
import meshi.util.file.MeshiWriter;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by chen on 14/09/2016.
 */
public class PostMortem12 {
    private static String target = "T0860";
    private static String rootDir = "D:\\Work\\CASP12\\postMortem\\";
    private static String outputFile = rootDir+target+"\\stage2.dat";
    private static String nativeFileName = rootDir+target+"\\"+target+".pdb";
    private static String scoreFileName = rootDir+target+"\\stage2\\ScoreFRes\\RESULTS."+target+".2.MESHI.SERVER.txt";
    private static String decoys = rootDir+target+"\\stage2\\pdb";
    public static void main(String[] args) throws IOException{
        Utils.verboseOff();
        Protein nativeStructure = new Protein(new AtomList(nativeFileName), new ResidueExtendedAtomsCreator());
        MeshiLineReader scoreReader = new MeshiLineReader(scoreFileName);
        MeshiWriter writer = new MeshiWriter(outputFile);
        String line;
        while ((line = scoreReader.readLine()) != null) {
            String[] words = line.split(" ");
            String decoyName = words[0];
            String decoyFileName = decoys + "\\" + decoyName + ".pdb";
            Protein decoy = new Protein(new AtomList(decoyFileName), new ResidueExtendedAtomsCreator());
            System.out.println(line + " " + Rms.gdt(nativeStructure, decoy)[0]);
            writer.println(line + " " + Rms.gdt(nativeStructure, decoy)[0]);
        }
        scoreReader.close();
        writer.close();
        int[] zvl = new int[7];
    }
}
