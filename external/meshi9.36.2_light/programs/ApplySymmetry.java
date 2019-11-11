package programs;

import meshi.geometry.Symmetry;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 30/05/2010...,,,,,,,,
 * Time: 18:59:42
 * To change this template use File | Settings | File Templates.
 */
public class ApplySymmetry {
    private static Protein protein;
    private static ArrayList<Symmetry> symmetryList;
    private static MeshiWriter output;

    public static void main(String[] argv) throws IOException {
        init(argv);


    }

    private static void init(String[] argv) throws IOException {
        if (argv.length != 2) {
            throw new RuntimeException(usageString(" Wrong number of arguments."));
        }
        String proteinFileName = argv[0];
        String outputFileName = proteinFileName.substring(0, proteinFileName.indexOf("pdb")) + "symmetry.pdb";
        output = new MeshiWriter(outputFileName);
        if (!(new File(argv[0])).exists()) throw new RuntimeException(usageString("Cannot find protein file."));
        AtomList tempAtomList = new AtomList(proteinFileName);
        new MolecularSystem();
        Protein protein = new Protein(tempAtomList, ResidueExtendedAtomsCreator.creator);

        String symmetryFileName = argv[1];
        if (!(new File(symmetryFileName).exists()))
            throw new RuntimeException(usageString("Cannot find protein file."));
        MeshiLineReader reader = new MeshiLineReader(symmetryFileName);
        symmetryList = Symmetry.getSymmetryList(reader);
    }

    private static String usageString(String message) {
        return "Error\n" + message + "\n" + "Usage: java -Xmx1000m programs/applySymmetry <protein file name> <Symmetry file>";
    }
}
