package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.sequences.AlignmentColumn;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.MeshiAttribute;

/**
 * Created by chen on 31/08/2017.
 */
public class CompareSequences {
    public static void main(String[] args)  throws AlignmentException{
        Protein protein1 = new Protein(new AtomList(args[0]), ResidueExtendedAtoms.creator);
        Protein protein2 = new Protein(new AtomList(args[1]), ResidueExtendedAtoms.creator);
        System.out.println("---------------------------------------\n"+
                            protein1.name()+" "+protein2.name());
        ResidueAlignment residueAlignment = new ResidueAlignment(protein1.chain(),"1",protein2.chain(),"2", ResidueAlignmentMethod.IDENTITY);
        System.out.println(residueAlignment.getSequenceAlignment()+"\n"+residueAlignment);
        for (AlignmentColumn column : residueAlignment) {
            Residue residue1 = (Residue) column.cell0().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            Residue residue2 = (Residue) column.cell1().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if (residue1.type() != residue2.type())
                System.out.println(residue1+" "+residue2);
        }
    }
}
