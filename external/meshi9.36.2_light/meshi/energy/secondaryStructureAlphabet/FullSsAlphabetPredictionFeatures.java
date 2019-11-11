package meshi.energy.secondaryStructureAlphabet;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.LocalStructure;
import meshi.sequences.MeshiSequence;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentColumn;
import meshi.util.MeshiAttribute;
import meshi.util.info.InfoType;

/**
 * Created by chen on 23/03/2016.
 */
public class FullSsAlphabetPredictionFeatures extends AbstractEnergy {


    private double[] structureAppearanceInNative = new double[LocalStructure.DsspSTRUCTURESortedByFrequency.length];
    private double[] correctPredictionCount = new double[LocalStructure.DsspSTRUCTURESortedByFrequency.length];
    private SequenceAlignment ssAlignment;

    MeshiSequence predictionSsSequence;
    MeshiSequence nativeSsSequence;


    public FullSsAlphabetPredictionFeatures(MeshiSequence predictionSsSequence, MeshiSequence nativeSsSequence) {
        super(toArray(), new StaticFeaturesInfo(), EnergyType.NON_DIFFERENTIAL);

        this.predictionSsSequence = predictionSsSequence;
        this.nativeSsSequence = nativeSsSequence;
        ssAlignment = new SequenceAlignment();
        ssAlignment.addAll(SequenceAlignment.identityAlignment(predictionSsSequence, nativeSsSequence));

        evaluate();
    }
    public EnergyInfoElement evaluate() {
        info.setValue(0.0);

        String str="param";
        for (int iParam=0; iParam < LocalStructure.DsspSTRUCTURESortedByFrequency.length; iParam++){

            ((StaticFeaturesInfo) info).params[iParam].setValue(predictionSuccessDSSP3(ssAlignment,new LocalStructure(LocalStructure.DsspSTRUCTURESortedByFrequency[iParam])));
        }
        return info;
    }

    private static double predictionSuccessDSSP3(SequenceAlignment alignment, LocalStructure dls){
        double sumPredicted = 0;
        double sumObserved  = 0;
        for (SequenceAlignmentColumn column : alignment){
            ResidueSsPrediction             predictionResidue            = (ResidueSsPrediction)             column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            LocalStructure nativeSsResidue         = (LocalStructure) column.cell1().getAttribute(MeshiAttribute.SS_DSSP);


            if (predictionResidue != null && nativeSsResidue != null && dls.compareTo(nativeSsResidue)==0) {
                sumObserved += predictionResidue.getProbabilityOf(nativeSsResidue.getDSSP3Reduction());
                sumPredicted++;
            }
        }
        if (sumPredicted == 0) return -1;
        else return sumObserved/sumPredicted;
    }




    private static double test1(SequenceAlignment alignment){
        return 0;
    }



    public void evaluateAtoms(){}
    public void test(TotalEnergy energy, Atom atom) {}





    private static class StaticFeaturesInfo extends EnergyInfoElement{
        public EnergyInfoElement[] params = new EnergyInfoElement[LocalStructure.DsspSTRUCTURESortedByFrequency.length];
        public StaticFeaturesInfo() {
            super(InfoType.SsS_PARAM, "secondaryStructureCompatibility");

            for (int iParam=0; iParam < params.length; iParam++){
                params[iParam]          = new EnergyInfoElement(InfoType.getTypeByTag(LocalStructure.DsspSTRUCTURESortedByFrequency[iParam]), "...");
                getChildren().add(params[iParam]);
            }

        }
    }

}
