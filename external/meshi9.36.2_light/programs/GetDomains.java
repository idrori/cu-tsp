package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ExtendedAtomsProtein;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.MeshiPotential;
import meshi.sequences.AlignmentColumn;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.UpdateableException;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;
import meshi.util.string.StringList;

import java.io.File;
import java.io.IOException;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 01/12/11
 * Time: 13:26
 * To change this template use File | Settings | File Templates.
 */
public class GetDomains extends MeshiProgram {
    public static void main(String[] argv) throws IOException, UpdateableException, AlignmentException{
        if (argv.length == 1)  {
            prepareCASP(argv[0]);
            return;
        }
        Utils.verboseOff();
        initRandom(0);
        if (argv.length > 0) getDomainsFromDirectories(argv);
    }

    public static void prepareCASP(String target) throws IOException{
        File serverModels = new File(target);
        serverModels.renameTo(new File(target+".serverModels"));
        serverModels = new File(target+".serverModels");
        File targetDir1 = new File(target);
        targetDir1.mkdir();
        File targetDir2 = new File(targetDir1.getName()+"\\"+target);
        targetDir2.mkdir();
        String newDir = targetDir1.getName()+"\\"+targetDir2.getName()+"\\serverModels";
        System.out.println(newDir);
        serverModels.renameTo(new File(newDir));
        File domainsDir = new File(targetDir1.getName()+"\\domains");
        domainsDir.mkdir();
        serverModels = new File(newDir);
        File[] files = serverModels.listFiles();
        long maxLength = 0;
        File maxFile = null;
        for (File file : files){
            if (file.length() > maxLength) {
                maxLength = file.length();
                maxFile = file;
            }
        }
        Files.copy(maxFile.toPath(), (new File(domainsDir.getAbsolutePath()+"\\"+target+"-D0.N.pdb")).toPath());
        files[0].renameTo(new File(domainsDir.getName() + "\\" + target + "-D0.N.pdb"));
    }

    public static void getDomainsFromDirectories(String[] argv) throws IOException, UpdateableException, AlignmentException{
        String prefix = argv[0];
        ArrayList <File>directories  = Utils.getDirectories(prefix);
        File domainsDir = getDomainsDir();
        File[] domainsFiles = domainsDir.listFiles();
        for (File file : domainsFiles) {
            String fileName =   file.getName();
            System.out.println(fileName);

            if(fileName.startsWith(prefix)) {
                for (File dir: directories) {
                    if (fileName.startsWith(dir.getName()))        {
                        MeshiWriter logWriter = new MeshiWriter(fileName+".log");
                        dir = new File(dir.getAbsolutePath()+"/serverModels");
                        if (dir.exists()) {
                            extractDomains(dir, file,logWriter,argv[1],argv[2].equals("filter")) ;
                /*            File newCopy = new File(dir, file.getName().substring(0,9)+"N.pdb");
                            StringList list = new StringList(file) ;
                            MeshiWriter writer = new MeshiWriter(newCopy);
                            list.print(writer);
                            writer.close();*/
                        }
                        else {
                            logWriter.println(dir.getAbsolutePath()+" does not exist.");
                        }
                        logWriter.close();
                    }
                }
            }
        }
    }

    public static void extractDomains(File dir, File domainFile,MeshiWriter logWriter, String commandsFile,Boolean filterFlag) throws IOException, AlignmentException, UpdateableException{
        Protein model;
        logWriter.println(dir.getAbsolutePath());
        System.out.println(dir.getAbsolutePath());
//        if (discontinuous(domainFile)) {
//            logWriter.println(domainFile+" is discontineous - ignored");
//            return;
//        }
        new MolecularSystem();
        Protein domain = ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(domainFile);
        new MolecularSystem();
        File domainDir = new File(domainFile.getName().substring(0,8));
        domainDir.mkdir();
        File[] files = dir.listFiles();
        for (File file:  files) {
            String newFileName;
            if ((model = getModel(domain, file, logWriter,commandsFile,filterFlag)) == null) continue;
            if (model.name().startsWith(domainFile.getName().substring(0, 9)))
                newFileName = domainFile.getName().substring(0, 9);
            else newFileName = domainFile.getName().substring(0, 9)+model.name();
            if (!newFileName.endsWith(".pdb"))
                newFileName +=".pdb";
            File newFile = new File(domainDir,newFileName);
            MeshiWriter outFile = new MeshiWriter(newFile);
            model.atoms().somewhere().print(outFile);
            outFile.close();
        }
        /*String newFileName = domainFile.getName().substring(0, 9)+"N.pdb";
        File newFile = new File(domainDir,newFileName);
        MeshiWriter outFile = new MeshiWriter(newFile);
        domain.getAtoms().print(outFile);
        outFile.close();*/
    }

    public static Protein getModel(Protein domain, File modelFile, MeshiWriter logWriter, String commandsFile,boolean filterFlage) throws UpdateableException, AlignmentException{

        ResidueList newResidueList = new ResidueList();
        double domainRadius = domain.atoms().radius();
        Protein model;

        new MolecularSystem();

        try {
            model =  ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(modelFile);
            logWriter.println(model.name());
            System.out.println(model.name());
        }
        catch (Exception ex) {
            return null;
        }
        if (numberOfGaps(model)>0) {
            logWriter.println("Model "+modelFile.getName()+" ignored - "+numberOfGaps(model)+" gaps");
            return null;
        }
        ResidueAlignment alignment;
        try {
            alignment = new ResidueAlignment(domain.chain(),domain.name(), model.chain(), model.name(), ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
        }
        catch(AlignmentException ex) {
            logWriter.println(ex);
            return null;
        }

        if (alignment.size() == 0) return null;

        int firstResidueNumber = ((Residue)alignment.get(0).cell1().object()).number();
        int lastResidueNumber =  ((Residue)alignment.get(alignment.size()-1).cell1().object()).number();
        //for (int iResidue = firstResidueNumber; iResidue <= lastResidueNumber; iResidue++)
        for (AlignmentColumn column : alignment) {
            newResidueList.add((Residue) column.cell1().object());
        }
        if (newResidueList.size()<= 10) {
            logWriter.println("Mode "+modelFile.getName()+" is too short "+newResidueList.size()+" ignored");
            return null;
        }
        newResidueList.setSourceFile(modelFile);
        AtomList newAtomList = newResidueList.atoms();
        double rg = newAtomList.radius();
        if ((rg > domainRadius*1.5) ||
                    ( 0.9*domain.chain().numberOfNonDummyResidues() > newResidueList.size())) {
            logWriter.print("Model "+modelFile.getName()+
                    " RG = "+rg+", Should be "+domainRadius+
                    " length "+newResidueList.size()+" should be "+domain.chain().numberOfNonDummyResidues());
            if (filterFlage) {
                    logWriter.println(" - ignored");
                    return null;
            }
            else{
                    logWriter.println("");
            }
        }

        Protein newModel = new Protein(newAtomList, new ResidueExtendedAtomsCreator());
        CommandList commands = new CommandList(commandsFile);

//        try {
//            Utils.addAtoms(newModel, true, commands,
//                    new PlaneParametersList(EnergyCreator.parametersDirectory(commands) +
//                            "/" + MeshiPotential.PLANE_PARAMETERS),false);
////        }
//        catch (Exception ex) {
//            logWriter.println("Exception was caught "+ex.getMessage()+" ignored");
//        }
        newModel.setName(model.name());
        return newModel;
    }


    public static boolean discontinuous(File proteinFile){
        new MolecularSystem();
        Protein protein = ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(proteinFile);
        return discontinuous(protein);
    }

    private static boolean discontinuous(Protein protein){
        ResidueList nonDummyResidues = protein.chain().getNonDummyResidues();
        int nNonDummyResidues = nonDummyResidues.size();
        if (nNonDummyResidues <0.9*(nonDummyResidues.get(nNonDummyResidues-1).number()-
                nonDummyResidues.get(0).number())) return true;
        return false;
    }
    public static int numberOfGaps(Protein model) {
        int i = 0;
        int numberOfGaps = 0;
        int nResidues = model.chain().size();
        while (model.chain().get(i).dummy())
            if (i < nResidues) i++;
            else throw new RuntimeException("weird");
        while (i < nResidues) {
            if (model.chain().get(i).dummy()) {
                numberOfGaps++;
                System.out.println("Gap in "+model+" in position "+i);
            }
            i++;
        }
        return numberOfGaps;
    }

    public static File getDomainsDir() {
        File thisDir = new File(".");
        File[] files = thisDir.listFiles();
        for (File file : files) {
            if (file.getName().equals("domains"))
                return file;
        }
        return null;
    }
}
