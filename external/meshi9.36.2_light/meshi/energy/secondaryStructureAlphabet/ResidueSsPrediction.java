package meshi.energy.secondaryStructureAlphabet;

import meshi.parameters.DsspLocalStructureLetter;
import meshi.parameters.ResidueType;
import meshi.util.MeshiAttribute;

/**
 * Created by chen on 23/03/2016.
 */
public class ResidueSsPrediction implements MeshiAttribute{
    ResidueType residueType;
    DsspLocalStructureLetter secondaryStructure;
    double helixProbability;
    double extendedProbability;
    double coilProbability;
    double highestProbability;

    public double getHighestProbability() {
        return highestProbability;
    }
    public ResidueType getResidueType() {return residueType;}

    public double getProbabilityOf(DsspLocalStructureLetter secondaryStructure) {
        if (secondaryStructure == DsspLocalStructureLetter.HELIX) return helixProbability;
        if (secondaryStructure == DsspLocalStructureLetter.SHEET) return extendedProbability;
        if (secondaryStructure == DsspLocalStructureLetter.COIL) return extendedProbability;
        return coilProbability;
    }

    public ResidueSsPrediction(String line) {
        String[] words;
        words = line.split("\\s+");

        residueType         = ResidueType.type(words[2]);
        if (words[3].charAt(0) == 'C') secondaryStructure = DsspLocalStructureLetter.COIL;
        else secondaryStructure  = DsspLocalStructureLetter.dsspLocalStructureLetter(1,words[3].charAt(0));
        coilProbability     = Double.valueOf(words[4]);
        helixProbability    = Double.valueOf(words[5]);
        extendedProbability = Double.valueOf(words[6]);
        highestProbability  = helixProbability;
        if (extendedProbability > highestProbability) highestProbability = extendedProbability;
        if (coilProbability     > highestProbability) highestProbability = coilProbability;
    }


    public int key() {
        return MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE;
    }




    public String toString() {
        String out = "ResidueSsPrediction "+residueType+" "+secondaryStructure+" "+helixProbability+" "+extendedProbability+" "+coilProbability;
        return out;

    }

}
