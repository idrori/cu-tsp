package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.MultiModelEnsemble;
import meshi.molecularElements.Protein;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.Rms;
import meshi.util.Utils;

import java.io.IOException;

/**
 * Created by chen on 12/09/2016.
 */
public class CalcGDT {
    public static void main(String[] args) throws IOException, AlignmentException{
        double gdt_ts;
        Utils.verboseOff();
        MultiModelEnsemble nativeEnsemble = new MultiModelEnsemble(args[0]);
        MultiModelEnsemble decoysEnsemble = new MultiModelEnsemble(args[1]);
        int modelNumber = 1;
        for (Protein decoy : decoysEnsemble) {
            gdt_ts = 0;
            String chainName = null;
            for (Protein nativeStructure : nativeEnsemble) {
                Object[] temp = calculateGDT(nativeStructure,decoy);
                gdt_ts += ((Double) temp[0]).doubleValue();
                chainName = (String) temp[1];
            }
            System.out.println(args[0]+" "+chainName+" "+args[1]+" model "+modelNumber+" "+gdt_ts/nativeEnsemble.size());
            modelNumber++;
        }
    }

    public static Object[] calculateGDT(Protein nativeStructure, Protein decoy) throws AlignmentException{
        ResidueAlignment bestAlignment = null;
        int chainLength = -10000;
        int longest = -100000;
        String chainName = "Weird";
        for (Chain chain :nativeStructure.chains()) {
            ResidueAlignment alignment = new ResidueAlignment(chain,"native",decoy.chain(),"decoy", ResidueAlignmentMethod.IDENTITY);
            if (alignment.size() > longest) {
                longest = alignment.size();
                bestAlignment = alignment;
                chainLength = chain.numberOfNonDummyResidues();
                chainName = chain.name();
            }
        }
        double[] gdt = Rms.gdt(bestAlignment,chainLength);
        //System.out.println(nativeStructure.name()+" "+decoy.name()+" "+gdt[0]);
        Object[] out = new Object[2];
        out[0] = new Double(gdt[0]);
        out[1] = chainName;
        return out;
    }
}
