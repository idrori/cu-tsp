package meshi.energy.secondaryStructureAlphabet;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.sequences.MeshiSequence;
import meshi.util.CommandList;
import meshi.util.dssp.DSSPFull;
import meshi.util.info.InfoType;

/**
 * Created by chen on 23/03/2016.
 */
public class SecondaryStructureAlphabetFeaturesCreator extends EnergyCreator{
    MeshiSequence[] SsSequences;

    public SecondaryStructureAlphabetFeaturesCreator() {
        super(InfoType.SECONDARY_STRUCTURE_COMPATIBILITY);


    }
    // The protein model is already assigned a DSSPFull structure for each residue
    public AbstractEnergy createEnergyTerm(Protein model, DistanceMatrix dm, CommandList commands) {
        String ssPredictionFileName = getNativeDsspFileName(commands);

        DSSPFull dssp = new DSSPFull(ssPredictionFileName);

        SsSequences = getDsspMeshiSequence(dssp); //takes the first chain sequence

        term =  new SecondaryStructureAlphabetFeatures(model, SsSequences);
        return term;

    }

    private static String getNativeDsspFileName(CommandList commands) {
        return commands.firstWord("nativeDsspFileName").secondWord();
    }

    private static MeshiSequence[] getDsspMeshiSequence(DSSPFull dssp){
        if (dssp.length() == 0) throw new RuntimeException("problem wth DSSP file");
        return dssp.getResidueSequenceWithFullSs();

    }


}
