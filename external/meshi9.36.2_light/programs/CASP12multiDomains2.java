package programs;

import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomType;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by chen on 25/05/2016.
 */
public class CASP12multiDomains2 {
    private static ArrayList<ResidueList> domains = new ArrayList<ResidueList>();
    private static ArrayList<String> fileNames = new ArrayList<String>();
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
            createModel(iModel, template, writer);
            //Utils.relax(template.atoms(), model, energyCreators1, commands);
            printModel(iModel, template, writer);
         }
        writer.close();

    }

    public static void printModel(int iModel, Protein model, MeshiWriter writer) {
        writer.println("MODEL " + iModel);
        writer.println("PARENT N/A");
        for (int iDomain = 0; iDomain < domains.size(); iDomain++) {
            ResidueList domain = domains.get(iDomain);
            writer.println("REMARK domain " + iDomain + " residues " + domain.get(0).number() + " - " + domain.get(domain.size() - 1).number() +
                    "extracted from " + fileNames.get(iDomain));
        }
        model.atoms().print(writer);
        writer.println("TER");
        writer.println("END");

    }

        public static void createModel(int iModel, Protein template, MeshiWriter writer){
        domains.clear();
        fileNames.clear();

        template.atoms().setNowhere();
        File modelDir = new File("model."+iModel+".domains");
        System.out.println("Extracting files from "+modelDir);
        for (String fileName : modelDir.list()) {
            if (fileName.endsWith("pdb")) {
                fileNames.add(fileName);
                System.out.println("    Added "+fileName);
            }
        }
        for (int iDomain = 1 ; iDomain <= fileNames.size(); iDomain++ ){
            for (String fileName : fileNames) {
                String domain = fileName.substring(6,8);
                if (domain.equals("D" + iDomain)) {
                    addDomain(modelDir.getAbsolutePath() + "/" + fileName, template, writer);
                }
            }
        }
        Atom prevDomainLastCA = null;
        for (ResidueList domain :domains) {
            if (prevDomainLastCA != null) {
                Atom firstCA = domain.get(0).ca();
                double dx = prevDomainLastCA.x() - firstCA.x() + 3.8;
                double dy = prevDomainLastCA.y() - firstCA.y();
                double dz = prevDomainLastCA.z() - firstCA.z();
                for (Residue residue : domain) {
                    for (Atom atom : residue.getAtoms())  {
                        atom.setXYZ(atom.x()+dx, atom.y()+dy, atom.z()+dz);
                    }
                }
            }
            prevDomainLastCA = domain.get(domain.size()-1).ca();
        }

    }

    private static void addDomain(String fileName, Protein template, MeshiWriter writer) {
        Protein model = new Protein(new AtomList(fileName), ResidueExtendedAtoms.creator);
        ResidueAlignment alignment;
        ResidueList domain = new ResidueList();
        domains.add(domain);
        try {
            alignment = new ResidueAlignment(template.chains(), "template", model.chains(), "model", ResidueAlignmentMethod.IDENTITY);
            //alignment.print(writer);
        } catch (AlignmentException ex) {throw new RuntimeException(ex);}
        for (ResidueAlignmentColumn column : alignment) {
            Residue templateResidue = column.residue0();
            if (!templateResidue.ca().nowhere()) continue;
            Residue modelResidue = column.residue1();
            if (templateResidue.type != modelResidue.type) {
                throw new RuntimeException("This is weird.");
            }
            domain.add(templateResidue);
            AtomList templateAtoms = templateResidue.getAtoms();
            AtomList modelAtoms    = modelResidue.getAtoms();
            for (Atom templateAtom : templateAtoms) {
                if (! templateAtom.nowhere()) continue;
                for (Atom modelAtom : modelAtoms) {
                    if (modelAtom.nowhere()) continue;
                    boolean found = false;
                    if (modelAtom.type() == templateAtom.type()) {
                        found = true;
                    }
                    else {
                        if ((modelAtom.type() == AtomType.TRC) && (templateAtom.name().equals("C")))
                            found = true;
                        else if ((modelAtom.type() == AtomType.TRO) && (templateAtom.name().equals("O")))
                            found = true;
                        else if ((modelAtom.type() == AtomType.TRN) && (templateAtom.name().equals("N")))
                            found = true;
                    }

                    if (found) {
                        templateAtom.setXYZ(modelAtom.x(), modelAtom.y(), modelAtom.z());
                        templateAtom.setTemperatureFactor(modelAtom.temperatureFactor());
                    }
                }
            }
        }

    }



    private static String  getHeader(String headerFileName) throws IOException{
        MeshiLineReader headerReader = new MeshiLineReader(headerFileName);
        String out = headerReader.readLine();
        String line;
        while ((line = headerReader.readLine()) != null) {
            out += "\n"+line;
        }
        return out;
    }
}
