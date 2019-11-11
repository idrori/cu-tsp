package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.CommandList;
import meshi.util.Rms;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by chen on 18/07/2016.
 */
public class FindModel {
    static Protein[] topModels = new Protein[5];
    static MeshiWriter writer;
    public static void main(String[] args) throws IOException{

        String targetFileName = args[0];
        //CommandList commands = new CommandList("commands");
        Protein target = new Protein(new AtomList(targetFileName),ResidueExtendedAtomsCreator.creator);
        File dir = new File("allServerModels");
        File[] files = dir.listFiles();
        for (int i = 0; i < 5; i++)
            files = getTopModel(target, files, i);


        String headerFileName = "../refinmentHeader.txt" ;
        String groupCode, scoreType;
        groupCode = "6553-5783-5054";
            scoreType = "MESHI";

    writer = new MeshiWriter(target.name()+".submission.pdb");
    String genericHeader = CASP12multiDomains.getHeader(headerFileName);
    String header = String.format(genericHeader,target.name().substring(0,5),groupCode,scoreType);
    writer.println(header);
    for (int iModel = 1; iModel <= 5; iModel++) {
        printModel(iModel);
    }
    writer.close();

}
    public static void printModel(int iModel){
        writer.println("MODEL " + iModel);
        writer.println("PARENT N/A");
        writer.println("REMARK "+ topModels[iModel-1].name());
        topModels[iModel-1].atoms().print(writer);
        writer.println("TER");
        writer.println("END");
    }

    private static File[] getTopModel(Protein target, File[] files,int index) {
        Protein best = null;
        double bestScore = -1;
        int bestI = -1;
        File[] out = new File[files.length - 1];
        for (int i = 0; i < files.length; i++) {
            Protein candidate = new Protein(new AtomList(files[i].getAbsolutePath()), ResidueExtendedAtomsCreator.creator);
            double[] gdt = Rms.gdt(target, candidate);
            double score = gdt[0];
            System.out.println(candidate.name()+" "+score);
            if (score > bestScore) {
                bestScore = score;
                best = candidate;
                bestI = i;
            }
        }
        System.out.println("*************\nBest is "+best.sourceFile()+" score: "+bestScore+"\n****************");
        topModels[index] = best;
        int j = 0;
        for (int i = 0; i < files.length; i++) {
            if (i != bestI) {
                out[j] = files[i];
                j++;
            }
        }
        return out;
    }

}
