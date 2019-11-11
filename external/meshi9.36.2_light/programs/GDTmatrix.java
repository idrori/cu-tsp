package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.MultiModelEnsemble;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.Pro;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.GDTcalculator;
import meshi.util.Rms;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by chen on 12/09/2016.
 */
public class GDTmatrix {
    private static String padding = "                                                                                ";
    public static void main(String[] args) throws IOException, AlignmentException{
        MultiModelEnsemble decoysEnsemble = new MultiModelEnsemble(args[0]);
        MeshiWriter writer = new MeshiWriter(args[1]);
        double[] gdt_tsTable = new double[decoysEnsemble.size()];
        double[] gdt_ts1Table = new double[decoysEnsemble.size()];
        double[] gdt_ts2Table = new double[decoysEnsemble.size()];
        double[] gdt_ts4Table = new double[decoysEnsemble.size()];
        double[] gdt_ts8Table = new double[decoysEnsemble.size()];
        double[] gdt_haTable = new double[decoysEnsemble.size()];
        double[] rmsTable = new double[decoysEnsemble.size()];
        double[] gdt_tsTable2 = new double[decoysEnsemble.size()];
        for (int i = 0; i < decoysEnsemble.size(); i++) {
            System.out.println(i+" of "+decoysEnsemble.size());
            Protein proteinI = decoysEnsemble.get(i);
            for (int j = i + 1; j < decoysEnsemble.size(); j++) {
                Protein proteinJ = decoysEnsemble.get(j);
                ResidueAlignment alignment = new ResidueAlignment(proteinI.chain(), "protein1",
                        proteinJ.chain(), "protein2", ResidueAlignmentMethod.IDENTITY);
                double rms = Rms.rms(alignment, Rms.RmsType.CA);
                double[] gdt1 = Rms.gdt(alignment, proteinI.chain().numberOfNonDummyResidues());
                double[] gdt2 = Rms.gdt(alignment, proteinJ.chain().numberOfNonDummyResidues());
                gdtMean(gdt_tsTable, i, j, gdt1[0], gdt2[0]);
                double avg = (gdt1[0]+gdt2[0])/2;
                gdtMean(gdt_tsTable2, i, j, avg*avg, avg*avg);
                gdtMean(gdt_ts1Table, i, j, gdt1[1], gdt2[1]);
                gdtMean(gdt_ts2Table, i, j, gdt1[2], gdt2[2]);
                gdtMean(gdt_ts4Table, i, j, gdt1[3], gdt2[3]);
                gdtMean(gdt_ts8Table, i, j, gdt1[4], gdt2[4]);
                gdt1 = Rms.gdt(alignment, proteinI.chain().numberOfNonDummyResidues(), GDTcalculator.Type.HA);
                gdt2 = Rms.gdt(alignment, proteinJ.chain().numberOfNonDummyResidues(), GDTcalculator.Type.HA);
                gdtMean(gdt_haTable, i, j, gdt1[0], gdt2[0]);
                rmsTable[i] += rms;
                rmsTable[j] += rms;
            }

        }
        Protein nativeStructure = null;
        if (args.length == 3) {
            nativeStructure = new Protein(new AtomList(args[0]),ResidueExtendedAtoms.creator);
        }
        int n = decoysEnsemble.size();
        for (int i = 0; i < decoysEnsemble.size(); i++) {
            String name = decoysEnsemble.get(i).name();
            String out = name+padding.substring(0,padding.length()-name.length())+" , "+
                    (gdt_tsTable[i]/n)+" , "+
                    (Math.sqrt(gdt_tsTable2[i]/n))+" , "+
                    (gdt_ts1Table[i]/n)+" , "+
                    (gdt_ts2Table[i]/n)+" , "+
                    (gdt_ts4Table[i]/n)+" , "+
                    (gdt_ts8Table[i]/n)+" , "+
                    (gdt_haTable[i]/n)+" , "+
                    (rmsTable[i]/n);
            if (nativeStructure != null) {
                ResidueAlignment alignment = new ResidueAlignment(nativeStructure.chain(), "nativeStructure",
                                        decoysEnsemble.get(i).chain(), "decoy", ResidueAlignmentMethod.IDENTITY);
                double[] gdt = Rms.gdt(alignment, nativeStructure.chain().numberOfNonDummyResidues());
                out += " , "+gdt[0];
            }
            System.out.println(out);
            writer.println(out);
        }
        writer.close();
    }

    public static void gdtMean(double[] table, int i, int j, double gdt1, double gdt2) {
        double mean = (gdt1 + gdt2)/2;
        table[i] += mean;
        table[j] += mean;
    }
}
