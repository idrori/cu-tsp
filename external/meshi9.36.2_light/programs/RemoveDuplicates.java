package programs;

import meshi.molecularElements.MultiModelEnsemble;
import meshi.molecularElements.Protein;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.GDTcalculator;
import meshi.util.Rms;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by chen on 12/09/2016.
 */
public class RemoveDuplicates {
    public static void main(String[] args) throws IOException, AlignmentException{
        MultiModelEnsemble decoysEnsemble = new MultiModelEnsemble(args[0]);
        for (int i = 0; i < decoysEnsemble.size(); i++) {
            System.out.println(i+" of "+decoysEnsemble.size());
            Protein proteinI = decoysEnsemble.get(i);
            if (proteinI == null)
                continue;
            for (int j = i + 1; j < decoysEnsemble.size(); j++) {
                Protein proteinJ = decoysEnsemble.get(j);
                if (proteinJ == null) continue;
                ResidueAlignment alignment = new ResidueAlignment(proteinI.chain(), "protein1",
                        proteinJ.chain(), "protein2", ResidueAlignmentMethod.IDENTITY);
                double rms = Rms.rms(alignment, Rms.RmsType.CA);
                if (rms < 0.00000001) {
                    System.out.println("Duplicates found "+proteinI+" "+proteinJ);
                    decoysEnsemble.set(j, null);
                }
            }
        }
        MeshiWriter meshiWriter = new MeshiWriter(args[1]);
        decoysEnsemble.print(meshiWriter);
        meshiWriter.close();
    }

}
