package meshi.energy.compatebility;

import meshi.parameters.ResidueType;
import meshi.parameters.SecondaryStructure;
import meshi.util.MeshiAttribute;
import meshi.util.Utils;

/**
 * Created by chen on 23/03/2016.
 */
public class ResidueSsPrediction implements MeshiAttribute{
    ResidueType residueType;
    SecondaryStructure secondaryStructure;
    double helixProbability;
    double extendedProbability;
    double coilProbability;
    double highestProbability;

    public double getHighestProbability() {
        return highestProbability;
    }
    public ResidueType getResidueType() {return residueType;}

    public double getProbabilityOf(SecondaryStructure secondaryStructure) {
        if (secondaryStructure == SecondaryStructure.HELIX) return helixProbability;
        if (secondaryStructure == SecondaryStructure.SHEET) return extendedProbability;
        return coilProbability;
    }

    public ResidueSsPrediction(String line) {
        String[] words;
        words = line.split("\\s+");

        residueType         = ResidueType.type(words[2]);
        secondaryStructure  = SecondaryStructure.secondaryStructure(words[3].charAt(0));
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
