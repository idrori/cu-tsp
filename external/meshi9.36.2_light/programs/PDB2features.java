package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.compatebility.StaticFeaturesCreator;
import meshi.energy.contacts.ContactsCreator;
import meshi.energy.goap.GoapCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZ5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZStd5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rapdf.RapdfCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.secondaryStructureAlphabet.SecondaryStructureAlphabetFeaturesCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZStdPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore.RamachandranCoreCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZStdRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.solvation.AtomEnvironmentCreator;
import meshi.energy.solvation.AtomEnvironmentCreatorSS;
import meshi.energy.solvation.SolvationCreator;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.molecularElements.MultiModelEnsemble;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.Pro;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.OptimizerException;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.GDTcalculator;
import meshi.util.dssp.DSSP;
import meshi.util.dssp.DSSPFull;
import meshi.util.file.MeshiWriter;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfoElement;
import meshi.util.info.ProteinInfoOLd;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static meshi.util.KeyWords.MINIMIZE;

/**
 * Created by chen on 12/09/2016.
 */
public class PDB2features extends MeshiProgram{
    private static String padding = "                                                                                ";
    private static long startTime;

    public static void main(String[] args) throws IOException, AlignmentException, UpdateableException,EvaluationException, OptimizerException {
        startTime = (new Date()).getTime();
        String ensembleFile = getArgument("decoysEnsemble", args);
        String nativeStructureName = getArgument("nativeStructure", args);
        String commandsFile = getArgument("commands", args);
        String outFile = getArgument("outFile", args);
        String dsspDirName = getArgument("dsspDir", args);
        String pdbDirName = getArgument("pdbDir", args);

        MultiModelEnsemble decoysEnsemble = new MultiModelEnsemble(ensembleFile);
        decoysEnsemble = removeDuplicates(decoysEnsemble);
        Protein nativeStructure = new Protein(new AtomList(nativeStructureName), ResidueExtendedAtoms.creator);
        CommandList commands = new CommandList(commandsFile);

        File dsspDir = new File(dsspDirName);
        File[] dsspFiles = dsspDir.listFiles();
        File pdbDir = new File(pdbDirName);
        File[] pdbFiles = pdbDir.listFiles();

        LinkedHashMap<InfoType, Double>[] features = initHashTables(decoysEnsemble.size());
        minimizeDecoys(decoysEnsemble, features, nativeStructure, dsspFiles, pdbFiles, commands);
        GDTmatrixFeatures(decoysEnsemble, features);
        writeFeatures(ensembleFile, decoysEnsemble,features, outFile);
    }

    public static void writeFeatures(String ensembleFile, MultiModelEnsemble decoysEnsemble, LinkedHashMap<InfoType, Double>[] features, String outFileName) throws IOException{
        MeshiWriter writer = new MeshiWriter(outFileName);
        writer.print("EnsembleFile, decoyName, ");
        Set<InfoType> keys = features[0].keySet();
        Iterator<InfoType> keyIter = keys.iterator();
        while (keyIter.hasNext()) {
            writer.print(keyIter.next().tag);
            if (keyIter.hasNext())
                writer.print(", ");
            else
                writer.println();
        }
        for (int i = 0; i < features.length; i++){
            LinkedHashMap<InfoType, Double> featuresLine = features[i];
            String decoyName = decoysEnsemble.get(i).name();
            decoyName = decoyName+padding.substring(0,padding.length()-name.length())+", ";

            writer.print(ensembleFile+", "+decoyName);
            keyIter = keys.iterator();
            while (keyIter.hasNext()) {
                writer.print(featuresLine.get(keyIter.next()));
                if (keyIter.hasNext())
                    writer.print(", ");
                else
                    writer.println();
            }
        }
        writer.close();
    }


    public static void minimizeDecoys(MultiModelEnsemble multiModelEnsemble, LinkedHashMap<InfoType, Double>[] features,
                                      Protein nativeStructure, File[] dsspFiles, File[] pdbFiles, CommandList commands)
            throws UpdateableException, EvaluationException, OptimizerException, AlignmentException, IOException{
        EnergyCreator[] energyCreators = initEnergyCreators();
        for (int i = 0; i < multiModelEnsemble.size(); i++) {
            Protein decoy = multiModelEnsemble.get(i);
            File dsspFile = getFile(dsspFiles,decoy.name());
            DSSPFull dssp = new DSSPFull(dsspFile.getAbsolutePath());
            Utils.AssignFullDSSP(decoy, dssp);
            File originalFile = getFile(pdbFiles, decoy.name());
            Protein originalStructure = new Protein(new AtomList(originalFile.getAbsolutePath()),ResidueExtendedAtoms.creator);
            TotalEnergy energy = new TotalEnergy(decoy, energyCreators, commands, "Energy");
//            SteepestDecent steepestDecent = new SteepestDecent(energy, 0.5, 1000, 200, 0.00000001, 0.5, 1.5);
//                steepestDecent.run();
            LBFGS lbfgs = Utils.getLBFGS(energy, commands, MINIMIZE);
            lbfgs.run();
            ModelAnalyzer analyzer = new ModelAnalyzer(decoy, nativeStructure, originalStructure, energy, ResidueAlignmentMethod.IDENTITY);
            ProteinInfoOLd proteinInfo = analyzer.analyze("");
            for (MeshiInfoElement element : proteinInfo) {
                features[i].put(element.type, element.doubleValue());
            }
//            features[i].put(InfoType.CHANGE, Rms.rms(decoy, originalStructure, ResidueAlignmentMethod.IDENTITY));
            long time = (new Date()).getTime();
            System.out.println(decoy.name()+" done "+((time-startTime)/60000)+" minutes");
        }
    }

    private static File getFile(File[] files, String name) {
        int start = name.indexOf("MODELT0")+5;
        int end = name.indexOf(".scwrl");
        name = name.substring(start, end);
        for (File file : files)
            if (file.getName().indexOf(name) != -1)
                return file;
        throw new RuntimeException(name+" not found");
    }

    public static EnergyCreator[] initEnergyCreators() {
        TetherCreator tetherAllCreator = new TetherCreator(InfoType.TETHER_ALL);
        TetherCreator residueTetherCreator = new TetherCreator(InfoType.TETHER);
        HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
        SolvateCreatorHBforMinimization solvateCreator = new SolvateCreatorHBforMinimization();
        SolvationCreator solvationCreator = new SolvationCreator();
        AtomEnvironmentCreator atomEnvironmentCreator = new AtomEnvironmentCreator(solvationCreator);
        AtomEnvironmentCreatorSS atomEnvironmentCreatorSS = new AtomEnvironmentCreatorSS(solvationCreator);
        RgCreator rgCreator = new RgCreator();
        RapdfCreator rapdfCreator = new RapdfCreator();

        CompositePropensityCreator propensityCreator = new CompositePropensityCreator();
        CooperativeZPropensityCreator cooperativeZPropensityCreator = new CooperativeZPropensityCreator(propensityCreator);
        CooperativeZStdPropensityCreator cooperativeZStdPropensityCreator = new CooperativeZStdPropensityCreator(cooperativeZPropensityCreator);

        AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
        CooperativeZ5typesSummaCreator cooperativeZSummaCreator = new CooperativeZ5typesSummaCreator(summaCreator);
        CooperativeZStd5typesSummaCreator cooperativeZStdSummaCreator = new CooperativeZStd5typesSummaCreator(cooperativeZSummaCreator);
        ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
        RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
        RamachandranCoreCreator ramachandranCoreCreator = new RamachandranCoreCreator(ramachCreator);
        CooperativeZRamachandranCreator cooperativeZRamachandranCreator = new CooperativeZRamachandranCreator(ramachCreator);
        CooperativeZStdRamachandranCreator cooperativeZStdRamachandranCreator = new CooperativeZStdRamachandranCreator(cooperativeZRamachandranCreator);
        RamachandranCreator ramachandranCreator = new RamachandranCreator();
        HydrogenBondsPairsCreator hydrogenBondsPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);
        FlatRamachCreator flatRamachCreator = new FlatRamachCreator();    //ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
        PlaneCreator planeCreator = new PlaneCreator();
        GoapCreator goapCreator = new GoapCreator();
        final ContactsCreator contactsCreator = new ContactsCreator();

        EnergyCreator[] energyCreators = {
                new BondCreator(),
                new AngleCreator(),
                planeCreator,
                new OutOfPlaneCreator(),
                ramachandranCreator,
                ramachCreator,
                cooperativeZRamachandranCreator,
                cooperativeZStdRamachandranCreator,
                ramachandranCoreCreator,
                propensityCreator,
                cooperativeZPropensityCreator,
                cooperativeZStdPropensityCreator,
                summaCreator,
                excludedVolCreator,
                cooperativeZSummaCreator,
                cooperativeZStdSummaCreator,
                solvationCreator,
                atomEnvironmentCreator,
                atomEnvironmentCreatorSS,
                solvateCreator,
                hydrogenBondsCreator,
                hydrogenBondsPairsCreator,
                new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
                new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
                contactsCreator,
                rgCreator,
                goapCreator,
                residueTetherCreator,
                flatRamachCreator,
                tetherAllCreator,
                rapdfCreator,
                new StaticFeaturesCreator(),
                new SecondaryStructureAlphabetFeaturesCreator()
        };
        return energyCreators;
    }


    public static LinkedHashMap<InfoType, Double>[] initHashTables(int n) {
        LinkedHashMap<InfoType, Double>[] features = new java.util.LinkedHashMap[n];
        for (int i = 0; i < features.length; i++)
            features[i] = new LinkedHashMap<>();
        return features;
    }

    public static void GDTmatrixFeatures(MultiModelEnsemble decoysEnsemble, LinkedHashMap<InfoType, Double>[] features) throws AlignmentException{
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
        int n = decoysEnsemble.size();
        for (int i = 0; i < decoysEnsemble.size(); i++) {
            features[i].put(InfoType.GDT_CONSENSUS, gdt_tsTable[i]/n);
            features[i].put(InfoType.GDT1_CONSENSUS, gdt_ts1Table[i]/n);
            features[i].put(InfoType.GDT2_CONSENSUS, gdt_ts2Table[i]/n);
            features[i].put(InfoType.GDT4_CONSENSUS, gdt_ts4Table[i]/n);
            features[i].put(InfoType.GDT8_CONSENSUS, gdt_ts8Table[i]/n);
            features[i].put(InfoType.GDT_HA_CONSENSUS, gdt_haTable[i]/n);
            features[i].put(InfoType.RMS_CONSENSUS, rmsTable[i]/n);
        }
    }

    public static void gdtMean(double[] table, int i, int j, double gdt1, double gdt2) {
        double mean = (gdt1 + gdt2)/2;
        table[i] += mean;
        table[j] += mean;
    }

    public static MultiModelEnsemble removeDuplicates(MultiModelEnsemble decoysEnsemble) throws AlignmentException{
            for (int i = 0; i < decoysEnsemble.size(); i++) {
                System.out.println(i+" of "+decoysEnsemble.size());
                Protein proteinI = decoysEnsemble.get(i);
                if (proteinI == null)
                    continue;
                for (int j = i + 1; j < decoysEnsemble.size(); j++) {
                    Protein proteinJ = decoysEnsemble.get(j);
                    if (proteinJ == null) continue;
                    ResidueAlignment alignment = new ResidueAlignment(proteinI.chain(), "protein1",
                            proteinJ.chain(), "protein2", ResidueAlignmentMethod.IDENTITY);
                    double rms = Rms.rms(alignment, Rms.RmsType.CA);
                    if (rms < 0.00000001) {
                        System.out.println("Duplicates found "+proteinI+" "+proteinJ);
                        decoysEnsemble.set(j, null);
                    }
                }
            }
            MultiModelEnsemble out = new MultiModelEnsemble();
            for (Protein decoy : decoysEnsemble)
                if (decoy != null)
                    out.add(decoy);
            return out;
    }
}
