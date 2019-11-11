package meshi.energy.compatebility;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.sequences.AlignmentColumn;
import meshi.sequences.MeshiSequence;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentCell;
import meshi.util.CharAttribute;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.InfoType;

import java.io.IOException;

/**
 * Created by chen on 23/03/2016.
 */
public class StaticFeaturesCreator extends EnergyCreator{
    MeshiSequence SsSequence;
    MeshiSequence accPredictionSequence;
    MeshiSequence accPrediction;

    public StaticFeaturesCreator() {
        super(InfoType.SECONDARY_STRUCTURE_COMPATIBILITY);


    }

    public AbstractEnergy createEnergyTerm(Protein model, DistanceMatrix dm, CommandList commands) {
        String ssPredictionFileName = getSsPredictionFileName(commands);
        String accesibilityPredictionFileName = getAccesibilityPredictionFileName(commands);

        SsSequence = new MeshiSequence("Sequence of SS prediction");
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
                SsSequence.add(cell);
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
        String scwrlFilesDirectory;
        if (commands.keyExists("scwrlFiles")) {
            scwrlFilesDirectory = commands.firstWord("scwrlFiles").secondWord();
        }
        else scwrlFilesDirectory = ".";
        String voromqaFilesDirectory;
        if (commands.keyExists("voromqaFiles")) {
            voromqaFilesDirectory = commands.firstWord("voromqaFiles").secondWord();
            term =  new StaticFeatures(model, SsSequence, accPredictionSequence, scwrlFilesDirectory, voromqaFilesDirectory);
        }
        else term =  new StaticFeatures(model, SsSequence, accPredictionSequence, scwrlFilesDirectory, null);
        return term;

    }

    private static String getSsPredictionFileName(CommandList commands) {
        return commands.firstWord("predictedSsFileName").secondWord();
    }

    private static String getAccesibilityPredictionFileName(CommandList commands) {
        return commands.firstWord("predictedAccFileName").secondWord();
    }

}
