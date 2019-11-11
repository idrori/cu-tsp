package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.extendedAtoms.ExtendedAtomsProtein;
import meshi.sequences.AlignmentException;
import meshi.util.UpdateableException;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 29/08/13
 * Time: 12:43
 * To change this template use File | Settings | File Templates.
 */
public class RemoveIllegalModels {
    public static void main(String[] args) throws IOException,UpdateableException, AlignmentException {
        File here = new File(".");
        File[] files = here.listFiles();
        for (File dir : files) {
            if (dir.isDirectory()){
                Protein model,nativeStructure = null;
                MeshiWriter logWriter = new MeshiWriter(dir.getAbsoluteFile()+".log");
                System.out.println(dir.getAbsoluteFile());
                File[] proteinFiles = dir.listFiles();
                for (File proteinFile : proteinFiles) {
                    if (proteinFile.getName().endsWith(".N.pdb")){
                        nativeStructure =  ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(proteinFile);
                    }
                }
                if (nativeStructure == null)
                    throw new RuntimeException("This is weird");

                for (File proteinFile : proteinFiles) {
                    if ((!proteinFile.getName().endsWith(".N.pdb")) && (proteinFile.getName().endsWith(".pdb"))){
                        model =  GetDomains.getModel(nativeStructure, proteinFile, logWriter,args[0],true);
                        if (model == null)
                            proteinFile.renameTo(new File(dir.getAbsoluteFile()+"/illegal/"+proteinFile.getName()));
                    }
                }
                logWriter.close();
            }
        }
    }
}
