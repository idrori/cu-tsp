package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomType;
import meshi.sequences.*;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by chen on 25/05/2016.
 */
public class CASP12multiDomains {
    public static void main(String[] args) throws IOException{
        String target = args[0];
        String headerFileName = "../weFoldHeader.txt" ;
        File workingDirectory = new File(".");
        String groupCode, scoreType;

        if (workingDirectory.getAbsolutePath().contains("Shokoufeh")) {
            groupCode = "2567-0767-7367";
            scoreType = "MSC";
        }
        else {
            groupCode = "8988-5959-0269";
            scoreType = "MESHI";
        }

        MeshiWriter writer = new MeshiWriter(target+".multiDomain.submission."+scoreType+".txt");
        String genericHeader = getHeader(headerFileName);
        String header = String.format(genericHeader,target,groupCode,scoreType);
        writer.println(header);
        Protein template = new Protein(new AtomList(target+".template.pdb"), ResidueExtendedAtoms.creator);
        for (int iModel = 1; iModel <= 5; iModel++) {
            printModel(iModel,template,writer);
         }
        writer.close();

    }

    public static void printModel(int iModel, Protein template, MeshiWriter writer){
        writer.println("MODEL " + iModel);

        template.atoms().setNowhere();
        File modelDir = new File("model."+iModel+".domains");
        System.out.println("Extracting files from "+modelDir);
        ArrayList<String> fileNames = new ArrayList();
        System.out.println(modelDir);
        for (String fileName : modelDir.list()) {
            if (fileName.endsWith("pdb")) {
                fileNames.add(fileName);
                System.out.println("    Added "+fileName);
            }
        }
        for (int iDomain = 1 ; iDomain <= fileNames.size(); iDomain++ ){
            for (String fileName : fileNames) {
                String domain = fileName.substring(6,8);
                if (domain.equals("D"+iDomain)) {
                    writer.println("REMARK domain extracted from "+fileName);
                    writer.println("PARENT N/A");
                    printDomain(modelDir.getAbsolutePath() + "/" + fileName, template, writer);
                    writer.println("TER");
                }
            }
        }
        writer.println("END");

    }

    private static void printDomain(String fileName, Protein template, MeshiWriter writer) {
        Protein model = new Protein(new AtomList(fileName), ResidueExtendedAtoms.creator);
        ResidueAlignment alignment;
        try {
            alignment = new ResidueAlignment(template.chain(), "template", model.chain(), "model", ResidueAlignmentMethod.IDENTITY);
            //alignment.print(writer);
        } catch (AlignmentException ex) {throw new RuntimeException(ex);}
        for (ResidueAlignmentColumn column : alignment) {
            Residue templateResidue = column.residue0();
            if (!templateResidue.ca().nowhere()) continue;
            Residue modelResidue = column.residue1();
            if (templateResidue.type != modelResidue.type) {
                throw new RuntimeException("This is weird.");
            }
            AtomList templateAtoms = templateResidue.getAtoms();
            AtomList modelAtoms    = modelResidue.getAtoms();
            for (Atom templateAtom : templateAtoms) {
                if (! templateAtom.nowhere()) continue;
                for (Atom modelAtom : modelAtoms) {
                    if (modelAtom.nowhere()) continue;
                    if (modelAtom.type() == templateAtom.type()) {
                        templateAtom.setXYZ(modelAtom.x(), modelAtom.y(), modelAtom.z());
                        templateAtom.setTemperatureFactor(modelAtom.temperatureFactor());
                    }
                }
            }
            templateAtoms.print(writer);
        }

    }



    public static String  getHeader(String headerFileName) throws IOException{
        MeshiLineReader headerReader = new MeshiLineReader(headerFileName);
        String out = headerReader.readLine();
        String line;
        while ((line = headerReader.readLine()) != null) {
            out += "\n"+line;
        }
        return out;
    }
}
