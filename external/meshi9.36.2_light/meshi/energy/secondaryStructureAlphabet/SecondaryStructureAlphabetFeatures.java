package meshi.energy.secondaryStructureAlphabet;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.DsspLocalStructureLetter;
import meshi.parameters.LocalStructure;
import meshi.parameters.LocalStructureAlphabetType;
import meshi.sequences.*;
import meshi.util.MeshiAttribute;
import meshi.util.Utils;
import meshi.util.filters.Filter;
import meshi.util.info.InfoType;

/**
 * Created by chen on 23/03/2016.
 */
public class SecondaryStructureAlphabetFeatures extends AbstractEnergy {

    private static final Filter NON_GAP = new NonGapFilter();
    private static final Filter SS = new SsFilter();
    private SequenceAlignment ssAlignment;
    Protein model;
    MeshiSequence[] sequenceWithPredictions;
    ResidueSequence residueSequence;


    public SecondaryStructureAlphabetFeatures(Protein model, MeshiSequence[] sequenceWithPredictions) {
        super(toArray(), new StaticFeaturesInfo(), EnergyType.NON_DIFFERENTIAL);
        if (! model.homoOligoMer())
            throw new RuntimeException("The current version can handle only monomers and homo-oligomers");

        this.model = model;
        this.sequenceWithPredictions = sequenceWithPredictions;
        ssAlignment = new SequenceAlignment();
        if (sequenceWithPredictions.length !=  model.chains().size())
            throw new RuntimeException("This is weird. sequenceWithPredictions.length & model.chains().size() = "+ sequenceWithPredictions.length +"  "+  model.chains().size());
        for (int chainI = 0; chainI < model.chains().size(); chainI++){
            Chain chain = model.chains().get(chainI);
            residueSequence = chain.sequence();
            MeshiSequence sequenceWithPrediction = sequenceWithPredictions[chainI];
            ssAlignment.addAll(SequenceAlignment.identityAlignment(sequenceWithPrediction, residueSequence));
        }
        evaluate();
    }

    public EnergyInfoElement evaluate() {
        info.setValue(ssCometability(ssAlignment,LocalStructureAlphabetType.DSSP3));
        ((StaticFeaturesInfo) info).dssp7.setValue(ssCometability(ssAlignment,LocalStructureAlphabetType.DSSP7));
        ((StaticFeaturesInfo) info).str2.setValue(ssCometability(ssAlignment,LocalStructureAlphabetType.STR2));
        ((StaticFeaturesInfo) info).dssp33.setValue(ssCometability(ssAlignment,LocalStructureAlphabetType.DSSP33));
        ((StaticFeaturesInfo) info).dssp103.setValue(ssCometability(ssAlignment,LocalStructureAlphabetType.DSSP103));
        //((StaticFeaturesInfo) info).coverage.setValue(getCoverage(ssAlignment, sequenceWithPrediction));
        //((StaticFeaturesInfo) info).coverage2.setValue(getStrictCoverage(ssAlignment, sequenceWithPrediction));
        //((StaticFeaturesInfo) info).modelNatoms.setValue(new Double(model.atoms().size()));
        //((StaticFeaturesInfo) info).length.setValue(new Double(model.residues().size()));
        //((StaticFeaturesInfo) info).one.setValue(new Double(1));
        return info;
    }

    private static double ssCometability(SequenceAlignment alignment, LocalStructureAlphabetType localStructureAlphabetType)   {
        double sumPredicted = 0;
        double sumObserved  = 0;
        for (SequenceAlignmentColumn column : alignment){
            LocalStructure prediction         = (LocalStructure) column.cell0().getAttribute(MeshiAttribute.SS_DSSP);
            Residue             residue            = (Residue)             column.cell1().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);

            if (residue != null && prediction != null) {
                if (residue.getLocalStructure().equals(prediction,localStructureAlphabetType)) sumObserved++;
                sumPredicted++;

            }
//            if (sumObserved != sumPredicted){
//                Utils.printDebug("yyyyyyyyyyyy","xxxxxxxxxxxxxxx ");
//            }
        }


        return sumObserved/sumPredicted;
    }

    private static double getStrictCoverage(SequenceAlignment alignment, MeshiSequence targetSequence) {
        int first = 0;
        int last = 0;
        boolean firstFound = false;
        Filter filter = SS;
        for (int i = 0; i < targetSequence.size(); i++) {
            if (filter.accept(targetSequence.get(i))) {
                if (!firstFound) {
                    first = i;
                    firstFound = true;
                }
                last = i;
            }
        }
        MeshiSequence newSequence = new MeshiSequence("Sequence without predicted N- and C-terminal coils");
        for (int i = first; i <= last; i++)
            newSequence.add(targetSequence.get(i));
        return getCoverage(alignment, newSequence);
    }

    private static double getCoverage(SequenceAlignment alignment, MeshiSequence targetSequence) {
        int nonGapColumns = 0;
        Filter filter = NON_GAP;
        for(AlignmentColumn column : alignment)
            if (filter.accept(column)) nonGapColumns++;
        double coverage = (1.0 * nonGapColumns)/targetSequence.size();
        if (coverage > 1) coverage = 1;
        return coverage;
    }

    private static class NonGapFilter implements Filter {
        public boolean accept(Object obj) {
            AlignmentColumn column = (AlignmentColumn) obj;
            return  (!column.cell0().gap()) && (!column.cell1().gap());
        }
    }

    private static class SsFilter implements Filter {
        public boolean accept(Object obj) {
            AlignmentColumn column = (AlignmentColumn) obj;
            ResidueSsPrediction ssp = (ResidueSsPrediction) column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            DsspLocalStructureLetter ss = ssp.secondaryStructure;
            if (ss == DsspLocalStructureLetter.HELIX) return true;
            if (ss == DsspLocalStructureLetter.SHEET) return true;
            return false;
        }
    }

    public void evaluateAtoms(){}
    public void test(TotalEnergy energy, Atom atom) {}





    private static class StaticFeaturesInfo extends EnergyInfoElement{
        public EnergyInfoElement dssp7,str2,dssp33,dssp103, coverage, coverage2,  modelNatoms, one, length;
        public StaticFeaturesInfo() {
            super(InfoType.SECONDARY_STRUCTURE_DSSP3_COMPATIBILITY, "secondaryStructureCompatibility");
            dssp7          = new EnergyInfoElement(InfoType.SECONDARY_STRUCTURE_DSSP7_COMPATIBILITY, "...");
            str2          = new EnergyInfoElement(InfoType.SECONDARY_STRUCTURE_STR2_COMPATIBILITY, "...");
            dssp33          = new EnergyInfoElement(InfoType.SECONDARY_STRUCTURE_DSSP33_COMPATIBILITY, "...");
            dssp103          = new EnergyInfoElement(InfoType.SECONDARY_STRUCTURE_DSSP103_COMPATIBILITY, "...");
            //coverage          = new EnergyInfoElement(InfoType.COVERAGE, "...");
            //coverage2         = new EnergyInfoElement(InfoType.SS_COVERAGE, "The fraction of the target protein, not including unstructured termini, which is covered by a model");
            //modelNatoms       = new EnergyInfoElement(InfoType.N_ATOMS, "The number of atoms in the model");
            //one               = new EnergyInfoElement(InfoType.ONE, "one");
            //length            = new EnergyInfoElement(InfoType.LENGTH, InfoType.LENGTH.tag);

            getChildren().add(dssp7);
            getChildren().add(str2);
            getChildren().add(dssp33);
            getChildren().add(dssp103);

          //  getChildren().add(coverage);
          //  getChildren().add(coverage2);
            //getChildren().add(modelNatoms);
            //getChildren().add(length);
            //getChildren().add(one);

        }
    }

}
