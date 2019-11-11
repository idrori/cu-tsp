package meshi.energy.compatebility;

import meshi.energy.*;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.SecondaryStructure;
import meshi.scoringFunctions.Sigmoid;
import meshi.sequences.*;
import meshi.util.CharAttribute;
import meshi.util.MeshiAttribute;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;
import meshi.util.info.ChainsInfo;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;

import java.io.File;

/**
 * Created by chen on 23/03/2016.
 */
public class StaticFeatures  extends AbstractEnergy implements EvaluatesResidues{
    private final double SS_FRACTION_SATURATION_MIDPOINT = 0.4;
    private final double SS_FRACTION_SATURATION_SLOPE = 100;
    private final double SS_SATURATION_MIDPOINT = 0.5;
    private final double SS_SATURATION_SLOPE = 10;
    private static final double SECONDARY_STRUCTURE_COMPATIBILITY_LIMIT = 0.8;
    private static final double SASA_COMPATIBILITY_LIMIT = 0.8;
    private static final Filter NON_GAP = new NonGapFilter();
    private static final Filter NON_GAP_SS = new NonGapSsFilter();
    private static final Filter SS = new SsFilter();
    private SequenceAlignment ssAlignment;
    private SequenceAlignment accAlignment;
    private MeshiSequence sequenceWithAccPrediction;
    Protein model;
    MeshiSequence sequenceWithPrediction;
    ResidueSequence residueSequence;
    String scwrlFilesDirectory;
    //String voromqaFileDirectory;

    public StaticFeatures(Protein model, MeshiSequence sequenceWithPrediction, MeshiSequence sequenceWithAccPrediction, String scwrlFilesDirectory, String voromqaFilesDirectory) {
        super(toArray(), new StaticFeaturesInfo(), EnergyType.NON_DIFFERENTIAL);
        if (! model.homoOligoMer())
            throw new RuntimeException("The current version can handle only monomers and homo-oligomers");
        this.sequenceWithAccPrediction = sequenceWithAccPrediction;
        this.model = model;
        this.sequenceWithPrediction = sequenceWithPrediction;
        ssAlignment = new SequenceAlignment();
        accAlignment = new SequenceAlignment();
        for (Chain chain : model.chains()) {
            residueSequence = chain.sequence();
            ssAlignment.addAll(SequenceAlignment.identityAlignment(sequenceWithPrediction, residueSequence));
            accAlignment.addAll(SequenceAlignment.identityAlignment(sequenceWithAccPrediction, residueSequence));
        }
        this.scwrlFilesDirectory = scwrlFilesDirectory;
        //this.voromqaFileDirectory = voromqaFilesDirectory;
        evaluate();
    }


//    private static double getSSCoverage(SequenceAlignment alignment, MeshiSequence targetSequence) {
//        int first = 0;
//        int last = 0;
//        boolean firstFound = false;
//        Filter filter = SS;
//        for (int i = 0; i < targetSequence.size(); i++) {
//            if (filter.accept(targetSequence.get(i))) {
//                if (!firstFound) {
//                    first = i;
//                    firstFound = true;
//                }
//                last = i;
//            }
//        }
//        MeshiSequence newSequence = new MeshiSequence("Sequence without predicted N- and C-terminal coils");
//        for (int i = first; i <= last; i++)
//            newSequence.add(targetSequence.get(i));
//        return getCoverage(alignment, newSequence);
//    }

    private static double getCoverage(SequenceAlignment alignment, MeshiSequence targetSequence, Filter filter1, Filter filter2) {
        int length = targetSequence.size();
        if (filter2 != null)
            for (SequenceAlignmentColumn column : targetSequence)
                if (!filter2.accept(column)) length--;
        int nonGapColumns = 0;
        for(AlignmentColumn column : alignment) {
            if (filter1.accept(column)) nonGapColumns++;
        }
        double coverage = (1.0 * nonGapColumns)/length;
        if (coverage > 1) coverage = 1;
        return coverage;
        /*
        int startTarget = 1000, endTarget = -1000, startModel = 1000, endModel = -1000;
        for (int iResidue = 0; iResidue < alignment.size(); iResidue++) {
            if ((startTarget == 1000) && filter.accept(alignment.get(iResidue).cell0())) startTarget = iResidue;
            if ((endTarget == -1000)  &&  filter.accept(alignment.get(alignment.size() - iResidue - 1).cell0()))
                endTarget = alignment.size() - iResidue - 1;
            if ((startModel == 1000)  && (!alignment.get(iResidue).cell1().gap())) startModel = iResidue;
            if ((endModel == -1000)   && (!alignment.get(alignment.size() - iResidue - 1).cell1().gap()))
                endModel = alignment.size() - iResidue - 1;
        }
        double out = (1.0*endModel-startModel)/(endTarget-startTarget);
        if (out > 1) out = 1;
        return out;                             */
    }

    private static class NonGapFilter implements Filter {
        public boolean accept(Object obj) {
            AlignmentColumn column = (AlignmentColumn) obj;
            return  (!column.cell0().gap()) && (!column.cell1().gap());
        }
    }
    private static class NonGapSsFilter implements Filter {
        public boolean accept(Object obj) {
            return SS.accept(obj) & NON_GAP.accept(obj);
        }
    }


    private static class SsFilter implements Filter {
        public boolean accept(Object obj) {
            AlignmentColumn column = (AlignmentColumn) obj;
            ResidueSsPrediction ssp = (ResidueSsPrediction) column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            if (ssp == null) {
                if (Utils.isStrict()) throw new RuntimeException("This is weird.");
                else return false;
            }
            SecondaryStructure ss = ssp.secondaryStructure;
            if (ss == SecondaryStructure.HELIX) return true;
            if (ss == SecondaryStructure.SHEET) return true;
            return false;
        }
    }

    private double sasaCompatibility(ChainsInfo chainsInfo) {
        double sum = 0;
        for (SequenceAlignmentColumn column : accAlignment){
            CharAttribute charAttribute = (CharAttribute) column.cell0().getAttribute(MeshiAttribute.CHAR);
            if (charAttribute == null) {
                String errorMessage = "A problem with the accAlignment\n"+accAlignment+"\n"+"in "+column;
                if (Utils.isStrict()) throw new RuntimeException(errorMessage);
                else continue;
            }
            char    accPrediction         = charAttribute.character;
            Residue residue               = (Residue)        column.cell1().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if (residue != null) {
                double relativeAccessinility = residue.getRelativeAccesibility();
                if (relativeAccessinility > 1) relativeAccessinility = 1;
                if (accPrediction == 'b') {
                    if (chainsInfo != null)
                        chainsInfo.get(residue.getChainNumber()).get(residue.number()).add(new DoubleInfoElement(InfoType.SASA_COMPATIBILITY, "sataCompatibility",1 - relativeAccessinility));
                    sum += 1 - relativeAccessinility;
                }
                else {
                    sum += relativeAccessinility;
                    if (chainsInfo != null)
                        chainsInfo.get(residue.getChainNumber()).get(residue.number()).add(new DoubleInfoElement(InfoType.SASA_COMPATIBILITY, "sataCompatibility",1 - relativeAccessinility));
                }
            }
        }
        return sum/accAlignment.size();
    }

    public double scwrlEnergy() {
        String scwrlFileName = scwrlFilesDirectory+"/"+model.sourceFile()+".scwrl";

        if (!(new File(scwrlFileName)).exists()) return 0;
        double energy = 0;
        String line;
        try {
            MeshiLineReader reader = new MeshiLineReader(scwrlFileName);
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("Total minimal energy")){
                    String[] words = line.split(" ");
                    energy = Double.valueOf(words[words.length-1]);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
        return energy;
    }

   /* public double voromqaEnergy() throws IOException{
        if (voromqaFileDirectory == null) return 0;
        File dir = new File(voromqaFileDirectory);
        if (! dir.isDirectory()) throw new RuntimeException("voromqaFileDirectory "+voromqaFileDirectory+" must be a directory");
        String[] temp = dir.list();
        if (temp.length != 1) throw new RuntimeException("voromqaFileDirectory "+voromqaFileDirectory+" must have a single file.");
        MeshiLineReader voromqaFile = new MeshiLineReader(temp[0]);
        String line = "";
        while ((line = voromqaFile.readLine()) != null)  {

        }
        String voromqaFileName = voromqaFileDirectory+"/"+model.sourceFile()+".voromqa";
        double energy = 0;
        try {
            MeshiLineReader reader = new MeshiLineReader(voromqaFileName);
            line = reader.readLine();
            String[] words = line.split(" ");
            energy = Double.valueOf(words[1]);
        }
        catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
        return energy;
    }                  */


    public void evaluateAtoms(){}
    public void test(TotalEnergy energy, Atom atom) {}

    public EnergyInfoElement evaluate() {
        evaluateResidues(null);
        return info;
    }

    public void evaluateResidues(ChainsInfo chainsInfo) {
        info.setValue(ssCompatibility(ssAlignment, chainsInfo));
        ((StaticFeaturesInfo) info).saturatedSsCompatibility.setValue(saturatedSsCometability(ssAlignment));
        ((StaticFeaturesInfo) info).helixCompatibility.setValue(ssCompatibility(ssAlignment, SecondaryStructure.HELIX));
        ((StaticFeaturesInfo) info).sheetCompatibility.setValue(ssCompatibility(ssAlignment, SecondaryStructure.SHEET));
        ((StaticFeaturesInfo) info).sasaCompatibility.setValue(sasaCompatibility(chainsInfo));
        ((StaticFeaturesInfo) info).saturatedSasaCompatibility.setValue(saturatedSasaCometability());
       ((StaticFeaturesInfo) info).coverage.setValue(getCoverage(ssAlignment, sequenceWithPrediction,NON_GAP, null));
        ((StaticFeaturesInfo) info).coverage2.setValue(getCoverage(ssAlignment, sequenceWithPrediction,NON_GAP_SS, SS));
        ((StaticFeaturesInfo) info).modelNatoms.setValue(new Double(model.atoms().size()));
        ((StaticFeaturesInfo) info).length.setValue(new Double(model.residues().size()));
        ((StaticFeaturesInfo) info).one.setValue(new Double(1));
        ((StaticFeaturesInfo) info).scwrlEnergy.setValue(scwrlEnergy());
        //((StaticFeaturesInfo) info).voromqaEnergy.setValue(voromqaEnergy());
        evaluateFractions();
    }

    public void evaluateFractions() {
        // sigmoid parameters - maps [0 1] to [0 1] with 0.4 mapped to 0.5
        int sum = 0;
        double sumSS = 0;
        double sumHelix = 0;
        double sumSheet = 0;
        for (Residue residue : model.residues()) {
            if (!residue.dummy()) {
                sum++;
                if (residue.getSecondaryStructure() != SecondaryStructure.COIL)
                    sumSS += 1;
                if (residue.getSecondaryStructure() == SecondaryStructure.HELIX)
                    sumHelix += 1;
                if (residue.getSecondaryStructure() == SecondaryStructure.SHEET)
                    sumSheet += 1;
            }
        }
        double ssFraction = sumSS/sum;
        Sigmoid sigmoid = new Sigmoid(SS_SATURATION_MIDPOINT, SS_SATURATION_SLOPE);
        double helixF;
        if (sumHelix == 0) helixF = 0.05; // for numerical stability, in case sumSheet is also zero.
        else helixF = sigmoid.sigmoid(sumHelix/(sumSheet+sumHelix))*0.9+0.05;
        double sheetF;
        if (sumSheet == 0) sheetF = 0.05;
        else sheetF = sigmoid.sigmoid(sumSheet/(sumSheet+sumHelix))*0.9+0.05;
        sigmoid = new Sigmoid(SS_FRACTION_SATURATION_MIDPOINT, SS_FRACTION_SATURATION_SLOPE);
        ((StaticFeaturesInfo) info).secondaryStructureFraction.setValue(ssFraction);
        ((StaticFeaturesInfo) info).saturatedSSfraction.setValue(sigmoid.sigmoid(ssFraction));
        ((StaticFeaturesInfo) info).helixFraction.setValue(helixF);
        ((StaticFeaturesInfo) info).sheetFraction.setValue(sheetF);
    }


    private static double saturatedSsCometability(SequenceAlignment alignment)   {
        double unSaturated = ssCompatibility(alignment, new ChainsInfo());
        if (unSaturated <= SECONDARY_STRUCTURE_COMPATIBILITY_LIMIT) return unSaturated/SECONDARY_STRUCTURE_COMPATIBILITY_LIMIT;
        return 1.0;
    }
    private double saturatedSasaCometability()   {
        double unSaturated = sasaCompatibility(null);
        if (unSaturated <= SASA_COMPATIBILITY_LIMIT) return unSaturated/SASA_COMPATIBILITY_LIMIT;
        return 1.0;
    }
    private static double ssCompatibility(SequenceAlignment alignment, ChainsInfo chainsInfo)   {
        Utils.println("Calculating ssCompatibility, alignment = \n"+alignment);
        if ((chainsInfo != null) && (chainsInfo.size() == 0))
            chainsInfo = null;
        double sumPredicted = 0;
        double sumObserved  = 0;
        for (SequenceAlignmentColumn column : alignment){
            ResidueSsPrediction prediction         = (ResidueSsPrediction) column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            Residue             residue            = (Residue)             column.cell1().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if (residue != null) {
                SecondaryStructure secondaryStructure = residue.getSecondaryStructure();
                    if (prediction != null) {
                        int chainNumber = residue.getChainNumber();
                        double observed  = prediction.getProbabilityOf(secondaryStructure);
                        double predicted = prediction.getHighestProbability();
                        if (chainsInfo != null)
                            chainsInfo.get(chainNumber).get(residue.number()).add(new DoubleInfoElement(InfoType.SECONDARY_STRUCTURE_COMPATIBILITY, "ssCompatibility",observed/predicted));
                        sumObserved +=  observed;
                        sumPredicted += predicted;
                    }
                    else {
                        String errorMessage =  "Something is wrong about the secondary structure prediction. Apparently it has holes in it.\n"+
                                "Residue is "+residue;
                        if (Utils.isStrict()) throw new RuntimeException(errorMessage);
                        else Utils.println(errorMessage);
                    }
            }
        }
        return sumObserved/sumPredicted;
    }

    private static double ssCompatibility(SequenceAlignment alignment, SecondaryStructure targetSS)   {
        double sumPredicted = 0;
        double sumObserved  = 0;
        for (SequenceAlignmentColumn column : alignment){
            ResidueSsPrediction prediction         = (ResidueSsPrediction) column.cell0().getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
            Residue             residue            = (Residue)             column.cell1().getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if ((residue != null) && (prediction != null)) {
                SecondaryStructure secondaryStructure = residue.getSecondaryStructure();
                if ((secondaryStructure == targetSS) && (prediction.secondaryStructure == targetSS)) {
                    sumObserved += prediction.getProbabilityOf(targetSS);
                    sumPredicted += prediction.getProbabilityOf(targetSS);
                }
                else if ((secondaryStructure != targetSS) && (prediction.secondaryStructure != targetSS)) {
                        sumObserved += 1 - prediction.getProbabilityOf(targetSS);
                        sumPredicted += 1 - prediction.getProbabilityOf(targetSS);
                }
                else if (secondaryStructure == targetSS) {
                        sumObserved += prediction.getProbabilityOf(targetSS);
                        sumPredicted += 1 - prediction.getProbabilityOf(targetSS);
                }
                else {
                        sumObserved += 1 - prediction.getProbabilityOf(targetSS);
                        sumPredicted += prediction.getProbabilityOf(targetSS);
                    }
                }

            }
        return sumObserved/sumPredicted;
    }

    private static class StaticFeaturesInfo extends EnergyInfoElement{
// with voro        public EnergyInfoElement secondaryStructureFraction, saturatedSSfraction, helixCompatibility, helixFraction, sheetFraction, sheetCompatibility, saturatedSsCompatibility, sasaCompatibility, saturatedSasaCompatibility, coverage, coverage2,  modelNatoms, one, length, scwrlEnergy, voromqaEnergy;
        public EnergyInfoElement secondaryStructureFraction, saturatedSSfraction, helixCompatibility, helixFraction, sheetFraction, sheetCompatibility, saturatedSsCompatibility, sasaCompatibility, saturatedSasaCompatibility, coverage, coverage2,  modelNatoms, one, length, scwrlEnergy;
        public StaticFeaturesInfo() {
            super(InfoType.SECONDARY_STRUCTURE_COMPATIBILITY, "secondaryStructureCompatibility");
            saturatedSsCompatibility = new EnergyInfoElement(InfoType.SATURATED_SECONDARY_STRUCTURE_COMPATIBILITY, "Compatibility with a prediction of secondary structure. Assuming 80% compatibility is as high as one may expect. ");
            helixCompatibility = new EnergyInfoElement(InfoType.HELIX_COMPATIBILITY," compatibility with helix prediction");
            sheetCompatibility = new EnergyInfoElement(InfoType.SHEET_COMPATIBILITY," compatibility with sheet prediction");
            sasaCompatibility = new EnergyInfoElement(InfoType.SASA_COMPATIBILITY, "Compatibility with prediction of Solvent accessible Surface Area ");
            saturatedSasaCompatibility = new EnergyInfoElement(InfoType.SATURATED_SASA_COMPATIBILITY, "Compatibility with prediction of Solvent accessible Surface Area. Assuming 80% compatibility is as high as one may expect. ");
            secondaryStructureFraction = new EnergyInfoElement(InfoType.SECONDARY_STRUCTURE_FRACTION,"The fraction of secondary structure residues.");
            saturatedSSfraction = new EnergyInfoElement(InfoType.SATURATED_SECONDARY_STRUCTURE_FRACTION, "A non-linear version of SECONDARY_STRUCTURE_FRACTION");coverage          = new EnergyInfoElement(InfoType.COVERAGE, "The fraction of the target protein which is covered by a model");
            helixFraction = new EnergyInfoElement(InfoType.HELIX_FRACTION, "Fraction of alpha helices among SS elements");
            sheetFraction = new EnergyInfoElement(InfoType.SHEET_FRACTION, "Fraction of beta sheet among SS elements");
            coverage2         = new EnergyInfoElement(InfoType.SS_COVERAGE, "The fraction of the target secondary structure elements of the protein, which is covered by a model.");
            modelNatoms       = new EnergyInfoElement(InfoType.N_ATOMS, "The number of atoms in the model");
            one               = new EnergyInfoElement(InfoType.ONE, "one");
            length            = new EnergyInfoElement(InfoType.LENGTH, InfoType.LENGTH.tag);
            scwrlEnergy       = new EnergyInfoElement(InfoType.SCWRL, InfoType.SCWRL.tag);
            //voromqaEnergy       = new EnergyInfoElement(InfoType.VOROMQA, InfoType.VOROMQA.tag);

            getChildren().add(saturatedSsCompatibility);
            getChildren().add(helixCompatibility);
            getChildren().add(sheetCompatibility);
            getChildren().add(sasaCompatibility);
            getChildren().add(saturatedSasaCompatibility);
            getChildren().add(secondaryStructureFraction);
            getChildren().add(saturatedSSfraction);
            getChildren().add(helixFraction);
            getChildren().add(sheetFraction);
            getChildren().add(coverage);
            getChildren().add(coverage2);
            getChildren().add(modelNatoms);
            getChildren().add(scwrlEnergy);
            //getChildren().add(voromqaEnergy);
            getChildren().add(length);
            getChildren().add(one);

        }
    }

    public boolean evaluatesResidues(){
        return true;
    }

}
