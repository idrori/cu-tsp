package meshi.util.statistics;

import meshi.molecularElements.Protein;
import meshi.molecularElements.ProteinMetaData;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.ProteinAnalyzer;
import meshi.util.Rms;
import meshi.util.info.MeshiInfo;
import meshi.util.info.ProteinInfo;

import java.io.File;
import java.util.ArrayList;

/**
 * Created by chen on 30/12/2015.
 */
public
class GdtRmsAnalyzer implements ProteinAnalyzer {
    boolean nativeStructureSearched = false;
    Protein nativeStructure;

    public ProteinInfo analyze(Protein protein) {
        if (!nativeStructureSearched)
            nativeStructure = getNativeStructure(protein);
        if (nativeStructure == null)
            return dummyGdtRms(protein);
        else return gdtRms(protein, nativeStructure);
    }

    private Protein getNativeStructure(Protein protein) {
        String directoryName = (String) protein.metaData().get(ProteinMetaData.MetaDataKey.DIRECTORY);
        String target = (String) protein.metaData().get(ProteinMetaData.MetaDataKey.TARGET);
        String nativeFileName = directoryName + "\\" + target + ".N.pdb";
        File nativeFile = new File(nativeFileName);

        if (!nativeFile.exists()) return null;

        return new Protein(new AtomList(nativeFileName), ResidueExtendedAtomsCreator.creator);
    }

    private ProteinInfo gdtRms(Protein protein, Protein nativeStructure) {
        try {
            double[] gdt = Rms.gdt(nativeStructure, protein);
            double rms = Rms.rms(nativeStructure, protein, ResidueAlignmentMethod.IDENTITY);
            return gdtRms(protein, gdt[0], rms);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private ProteinInfo dummyGdtRms(Protein protein) {
        return gdtRms(protein, -1, -1);
    }

    private ProteinInfo gdtRms(Protein protein, double gdt, double rms) {
        ArrayList<MeshiInfo> infoList = new ArrayList();
        infoList.add(new MeshiInfo(null, new Double(gdt), "Gdt"));
        infoList.add(new MeshiInfo(null, new Double(rms), "Rms"));
        ProteinInfo proteinInfo = new ProteinInfo(protein.metaData(), infoList, protein);
        return proteinInfo;
    }


}
