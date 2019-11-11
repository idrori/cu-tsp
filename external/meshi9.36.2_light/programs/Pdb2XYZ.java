package programs;

import com.sun.org.apache.regexp.internal.RE;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.ResidueType;
import meshi.util.file.MeshiWriter;
import meshi.util.formats.Format;

import java.io.IOException;

public class Pdb2XYZ {
    public static void main(String[] args) throws IOException{
        Protein protein = new Protein(new AtomList(args[0]), ResidueExtendedAtoms.creator);
        MeshiWriter writer = new MeshiWriter(protein.name().substring(0,protein.name().length()-4)+".xyz");
        Chain chain = protein.chain();
        writer.println(chain.numberOfNonDummyResidues());
        writer.println(protein.name().substring(0,protein.name().length()-4)+"  CB");
        for (Residue residue : chain)
            if (! residue.dummy()) {
                if (residue.type == ResidueType.GLY)
                    writer.println(residue.type().nameOneLetter() + "\t" + residue.ca().x() + "\t" + residue.ca().y() + "\t" + residue.ca().z());
                else
                    writer.println(residue.type().nameOneLetter() + "\t" + residue.cb().x() + "\t" + residue.cb().y() + "\t" + residue.cb().z());
            }

        writer.close();
    }
}
