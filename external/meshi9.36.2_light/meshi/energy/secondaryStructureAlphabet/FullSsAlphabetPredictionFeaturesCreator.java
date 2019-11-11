package meshi.energy.secondaryStructureAlphabet;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.sequences.AlignmentColumn;
import meshi.sequences.MeshiSequence;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentCell;
import meshi.util.CharAttribute;
import meshi.util.CommandList;
import meshi.util.dssp.DSSPFull;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.InfoType;

import java.io.IOException;

/**
 * Created by chen on 23/03/2016.
 */
public class FullSsAlphabetPredictionFeaturesCreator extends EnergyCreator{
    MeshiSequence nativeSsSequence;
    MeshiSequence predictionSsSequence;
    MeshiSequence accPredictionSequence;
    MeshiSequence accPrediction;

    public FullSsAlphabetPredictionFeaturesCreator() {
        super(InfoType.SECONDARY_STRUCTURE_COMPATIBILITY);


    }

    public AbstractEnergy createEnergyTerm(Protein model, DistanceMatrix dm, CommandList commands) {
        String ssPredictionFileName = getSsPredictionFileName(commands);
        String accesibilityPredictionFileName = getAccesibilityPredictionFileName(commands);

        predictionSsSequence = new MeshiSequence("Sequence of SS prediction");
        MeshiLineReader reader = new MeshiLineReader(ssPredictionFileName);
        String line;
        int i = 0;
        String line1, line2;
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                if (line.length() < 10) continue;
                ResidueSsPrediction residueSsPrediction = new ResidueSsPrediction(line);
                SequenceAlignmentCell cell = new SequenceAlignmentCell(residueSsPrediction.residueType.nameOneLetter().charAt(0), i);
                i++;
                predictionSsSequence.add(cell);
                cell.addAttribute(residueSsPrediction);
            }

            reader.close();
            reader = new MeshiLineReader(accesibilityPredictionFileName);
            reader.readLine(); //title
            line1 = reader.readLine();
            line2 = reader.readLine();
            reader.close();
        } catch (IOException ex){
            ex.printStackTrace();
            throw new RuntimeException(ex.getMessage());
        }
        accPredictionSequence = new MeshiSequence(line1,"The sequence from accessability prediction");
        accPrediction         = new MeshiSequence(line2,"The prediction");
        SequenceAlignment alignment         = new SequenceAlignment(accPredictionSequence,accPrediction);
        for (final AlignmentColumn column : alignment) {
            char c = ((Character) column.cell1().obj).charValue();
            column.cell0().addAttribute(new CharAttribute(c));
        }

        String ssNativeDsspFileName = commands.firstWord("nativeDsspFileName").secondWord();
        DSSPFull nativeDssp = new DSSPFull(ssNativeDsspFileName);
        nativeSsSequence = getDsspMeshiSequence(nativeDssp)[0]; //takes the first chain sequence


        term =  new FullSsAlphabetPredictionFeatures(predictionSsSequence, nativeSsSequence);
        return term;

    }

    private static String getSsPredictionFileName(CommandList commands) {
        return commands.firstWord("predictedSsFileName").secondWord();
    }

    private static String getAccesibilityPredictionFileName(CommandList commands) {
        return commands.firstWord("predictedAccFileName").secondWord();
    }
    private static MeshiSequence[] getDsspMeshiSequence(DSSPFull dssp){
        if (dssp.length() == 0) throw new RuntimeException("problem wth DSSP file");
        return dssp.getResidueSequenceWithFullSs();

    }
}
