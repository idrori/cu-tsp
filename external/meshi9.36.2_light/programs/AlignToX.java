package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by chen on 07/01/2017.
 */
public class AlignToX {
    public static void main(String[] args) throws IOException {
        String inFileName = args[0];
        String outFileName = args[1];

        Protein inProtein = new Protein(new AtomList(inFileName), ResidueExtendedAtomsCreator.creator);
        Utils.alignToX(inProtein);
        inProtein.printAtomsToFile(outFileName);
    }
}
