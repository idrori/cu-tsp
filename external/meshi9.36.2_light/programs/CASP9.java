package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 24/05/2010
 * Time: 13:00:51
 * To change this template use File | Settings | File Templates.
 */
public class CASP9 {
    public static void main(String[] argv) throws IOException {
        MeshiWriter writer;
        String parentString1, parentString2, parentString3, parentString4, parentString5;
        AtomList atoms1, atoms2, atoms3, atoms4, atoms5;
        Protein model1, model2, model3, model4, model5;

        MeshiLineReader inFile = new MeshiLineReader(argv[0]);
        String target = inFile.readLine();

        String modelFileName1 = inFile.readLine();
        String modelFileName2 = inFile.readLine();
        String modelFileName3 = inFile.readLine();
        String modelFileName4 = inFile.readLine();
        String modelFileName5 = inFile.readLine();

        parentString1 = getParentString(modelFileName1);
        parentString2 = getParentString(modelFileName2);
        parentString3 = getParentString(modelFileName3);
        parentString4 = getParentString(modelFileName4);
        parentString5 = getParentString(modelFileName5);

        atoms1 = new AtomList(modelFileName1);
        atoms2 = new AtomList(modelFileName2);
        atoms3 = new AtomList(modelFileName3);
        atoms4 = new AtomList(modelFileName4);
        atoms5 = new AtomList(modelFileName5);

        new MolecularSystem();
        model1 = new Protein(atoms1, ResidueExtendedAtomsCreator.creator);

        new MolecularSystem();
        model2 = new Protein(atoms2, ResidueExtendedAtomsCreator.creator);

        new MolecularSystem();
        model3 = new Protein(atoms3, ResidueExtendedAtomsCreator.creator);

        new MolecularSystem();
        model4 = new Protein(atoms4, ResidueExtendedAtomsCreator.creator);


        new MolecularSystem();
        model5 = new Protein(atoms5, ResidueExtendedAtomsCreator.creator);

        writer = new MeshiWriter(target + ".casp9");
        printHeaders(writer, target);

        print(model1, 1, writer, parentString1);
        print(model2, 2, writer, parentString2);
        print(model3, 3, writer, parentString3);
        print(model4, 4, writer, parentString4);
        print(model5, 5, writer, parentString5);

        writer.close();
    }


    private static void printHeaders(MeshiWriter writer, String target) {
        writer.println("PFRMAT TS");
        writer.println("TARGET " + target);
        writer.println("AUTHOR 1208-2919-8492");
        writer.println("METHOD  All server models were energy minimized. The highest ranking models");
        writer.println("METHOD  were further optimized using MCM");
    }

    private static void print(Protein model, int number, MeshiWriter writer, String parentString) {
        writer.println("MODEL " + number);

        if (parentString != null) {
            writer.println(parentString);
        } else writer.println("PARENT  N/A");

        String originalModel = model.name().substring(0, model.name().indexOf("TS") + 3);
        writer.println("REMARK  This model is an attempt to energy optimize the server model " + originalModel);
        writer.println("REMARK Original file " + model.name());
        model.atoms().print(writer);
        writer.println("TER");
        writer.println("END\n");
    }

    private static String getParentString(String fileName) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String out;
        while ((out = reader.readLine()) != null) {
            if (out.startsWith("PARENT"))
                return out;
        }
        return null;
    }
}
