package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.MultiModelEnsemble;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.Rms;
import meshi.util.Utils;

import java.io.IOException;

/**
 * Created by chen on 12/09/2016.
 */
public class GDTcalculator {
    public static void main(String[] args) throws IOException, AlignmentException{
        Utils.verboseOff();
        Protein nativeStructure = new Protein(new AtomList(args[0]), ResidueExtendedAtoms.creator);
        Protein model = new Protein(new AtomList(args[1]), ResidueExtendedAtoms.creator);
        double gdt_ts = calculateGDT(nativeStructure,model);
        System.out.println("GDT_TS of "+nativeStructure+" "+model+" = "+gdt_ts);
    }

    public static double calculateGDT(Protein nativeStructure, Protein model) throws AlignmentException{
        ResidueAlignment alignment = new ResidueAlignment(nativeStructure.chain(),"native",
                                                          model.chain(),"model", ResidueAlignmentMethod.IDENTITY);
        double[] gdt = Rms.gdt(alignment,nativeStructure.chain().numberOfNonDummyResidues());
        return gdt[0];
    }
}
