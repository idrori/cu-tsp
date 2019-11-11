package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.*;
import meshi.util.MeshiAttribute;
import meshi.util.Utils;

/**
 * Created by chen on 09/01/2017.
 */
public class ForOfir {
    public static void main(String[] args) {
        String modelFileName = args[0];
        String dsspFileName = args[1];
        String ssPredictionFileName = args[2];

        Protein model = new Protein(new AtomList(modelFileName), ResidueExtendedAtomsCreator.creator);
        Utils.AssignDSSP(model, dsspFileName);
        MeshiSequence ssSequense = Utils.getSSPrediction(ssPredictionFileName);
        SequenceAlignment ssAlignment = new SequenceAlignment();
        ResidueSequence residueSequence = model.chain().sequence();
        ssAlignment.addAll(SequenceAlignment.identityAlignment(ssSequense, residueSequence));
        for (SequenceAlignmentColumn column : ssAlignment) {
            AlignmentCell cell0 = column.cell0();
            AlignmentCell cell1 = column.cell1();
            System.out.println( cell0.getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE).getClass());
            System.out.println( cell1.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE).getClass());
            System.out.println( cell0+" "+cell0.getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE)+" ; "+column.cell1()+" "+cell1.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE));
        }
    }
}
