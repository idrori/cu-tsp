/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.PDB.PdbLineATOM;
import meshi.PDB.PdbLineIterator;
import meshi.PDB.PdbReader;
import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.conservationContacts.ConservationContactsCreator;
import meshi.energy.conservationContacts.ConservationContactsHrCreator;
import meshi.energy.contacts.ContactsCreator;
import meshi.energy.goap.GoapCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.one.OneCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZ5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZStd5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rapdf.RapdfCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZStdPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore.RamachandranCoreCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZStdRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateType;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.solvation.AtomEnvironmentCreator;
import meshi.energy.solvation.AtomEnvironmentCreatorSS;
import meshi.energy.solvation.SolvationCreator;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ProteinMetaData;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.*;
import meshi.parameters.MeshiPotential;
import meshi.parameters.SecondaryStructure;
import meshi.scoringFunctions.CombinedEnergyScore;
import meshi.scoringFunctions.GdtScore;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.info.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

//import meshi.energy.one.OneCreator;

public class Evaluate extends MeshiProgram implements KeyWords {

    public static final String NAME = "Evaluate";

    private static String inFileName, outFileName, datFile, dsspFile;
    private static double startGdt_ha;
    private static Protein model, originalModel;
    private static Boolean OK;
    private static CommandList commands;
    private static TotalEnergy energy;
    private static ArrayList<Score> scoreFunctions;
    private static TotalEnergy perturbationEnergy1, perturbationEnergy2, perturbationEnergy3, perturbationEnergy4;
    private static DistanceMatrix distanceMatrix;
    private static EvaluateLogger log;
    private static int seed;


    private static TetherCreator tetherAllCreator = new TetherCreator(InfoType.TETHER_ALL);
    private static TetherCreator residueTetherCreator = new TetherCreator(InfoType.TETHER);
    private static HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
    private static SolvateCreatorHBforMinimization solvateCreator = new SolvateCreatorHBforMinimization();
    private static SolvationCreator solvationCreator = new SolvationCreator();
    private static AtomEnvironmentCreator atomEnvironmentCreator = new AtomEnvironmentCreator(solvationCreator);
    private static AtomEnvironmentCreatorSS atomEnvironmentCreatorSS = new AtomEnvironmentCreatorSS(solvationCreator);
    private static RgCreator rgCreator = new RgCreator();
    //private static ContactsCreator contactsCreator = new ContactsCreator(InfoType.CONTACTS_ETC);
    private static RapdfCreator rapdfCreator = new RapdfCreator();

    private static CompositePropensityCreator propensityCreator = new CompositePropensityCreator();
    private static CooperativeZPropensityCreator cooperativeZPropensityCreator = new CooperativeZPropensityCreator(propensityCreator);
    private static CooperativeZStdPropensityCreator cooperativeZStdPropensityCreator = new CooperativeZStdPropensityCreator(cooperativeZPropensityCreator);

    private static AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
    private static CooperativeZ5typesSummaCreator cooperativeZSummaCreator = new CooperativeZ5typesSummaCreator(summaCreator);
    private static CooperativeZStd5typesSummaCreator cooperativeZStdSummaCreator = new CooperativeZStd5typesSummaCreator(cooperativeZSummaCreator);
    private static ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
    private static RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
    private static RamachandranCoreCreator ramachandranCoreCreator = new RamachandranCoreCreator(ramachCreator);
    private static CooperativeZRamachandranCreator cooperativeZRamachandranCreator = new CooperativeZRamachandranCreator(ramachCreator);
    private static CooperativeZStdRamachandranCreator cooperativeZStdRamachandranCreator = new CooperativeZStdRamachandranCreator(cooperativeZRamachandranCreator);
    private static RamachandranCreator ramachandranCreator = new RamachandranCreator();
    private static InflateCreator inflateCreator = new InflateCreator(InflateType.SIMPLE);
    private static InflateCreator inflatePerSegmentCreator = new InflateCreator(InflateType.PER_SEGMENT);
    private static InflateCreator inflateBySegmentCreator = new InflateCreator(InflateType.BY_SEGMENT);
    private static ConservationContactsCreator conservationContactsCreator8 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS8);
    private static ConservationContactsCreator conservationContactsCreator11 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS11);
    private static ConservationContactsCreator conservationContactsCreator15 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS15);
    private static ConservationContactsHrCreator conservationContactsHrCreator = new ConservationContactsHrCreator();
    private static HydrogenBondsPairsCreator hydrogenBondsPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);
    private static  FlatRamachCreator flatRamachCreator = new FlatRamachCreator();    //private static ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
    //private static OneCreator oneCreator = new OneCreator();
    private static PlaneCreator planeCreator = new PlaneCreator();
    private static GoapCreator goapCreator = new GoapCreator();
    private static final ContactsCreator  contactsCreator  = new ContactsCreator();

    private Score optimizationScore = null;
    private static EnergyCreator[] energyCreators = {
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
            conservationContactsCreator8,
            conservationContactsCreator11,
            conservationContactsCreator15,
            conservationContactsHrCreator,
            new OneCreator(),
		    rapdfCreator,
    };

    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException, AlignmentException {
        init(argv);
        PdbReader pdbReader = new PdbReader(inFileName);
        MeshiLineReader datReader = new MeshiLineReader(datFile);
        datReader.readLine();
        datReader.readLine();
        String line = pdbReader.readLine();
        int modelNumber = 1;
        while ((line = pdbReader.readLine()) != null) {
            PdbLineIterator pdbLineIterator = new PdbLineIterator(pdbReader);
            if (!pdbLineIterator.hasNext()) break;
            model = new Protein(new AtomList(pdbLineIterator,new PdbLineATOM()),ResidueExtendedAtoms.creator,commands);
            addAtoms(model,commands);
            String datLine = datReader.readLine();
            String[] words = datLine.split(" ");
            double gdt_ha = Double.parseDouble(words[7]);
            int trIndex = datLine.indexOf("TR");
            String target = datLine.substring(trIndex, trIndex+5);
            model.setTarget(target);
            Utils.alignToX(model);
            Utils.AssignDSSP(model, dsspFile);
            Utils.setSS(model, commands);
            scoreFunctions = GdtScore.getScoreFunctions(commands);
            energy = new TotalEnergy(model,energyCreators, commands, "My energy");
            ModelAnalyzer analyzer = new ModelAnalyzer(model, null, originalModel, energy, ResidueAlignmentMethod.IDENTITY);
            ProteinInfoOLd proteinInfoOLd = analyzer.analyze("DATA");
            MeshiInfoElement gdt_ha_element = new DoubleInfoElement(InfoType.GDT_HA, "gdt_ha",gdt_ha);
            proteinInfoOLd.add(gdt_ha_element);
            MeshiInfoElement start_gdt_ha_element = new DoubleInfoElement(InfoType.START_GDT_HA, "start_gdt_ha",startGdt_ha);
            proteinInfoOLd.add(start_gdt_ha_element);
            MeshiInfoElement length_element = new DoubleInfoElement(InfoType.LENGTH, "start_gdt_ha",model.chain().numberOfNonDummyResidues());
            proteinInfoOLd.add(length_element);
            outFileName = inFileName+"."+modelNumber+".xml";
            modelNumber++;
            MeshiInfoXMLwriter outWriter = new MeshiInfoXMLwriter(outFileName);
            print(proteinInfoOLd, outWriter,1, target);
            outWriter.close();
            /*
            ProteinInfoList proteinInfoList= new ProteinInfoList();
            ArrayList<MeshiInfo> infoList = new ArrayList();
            infoList.add(energy.energyInfo());
            ProteinInfo info = new ProteinInfo(model.metaData(),infoList, model);
            info.add(info);      */
            Utils.printDebug("main ",model.name() + datLine);
        }

            //((PlaneEnergy)planeCreator.term()).fixCisTrans();


    }
    public static void print(ProteinInfoOLd proteinInfo, MeshiInfoXMLwriter writer, int indentationTabs, String target) throws IOException {
        writer.println("<ProteinInfoList name=\"" + target + "\">");
        writer.println(proteinInfo.toXML(indentationTabs + 1));
        writer.println("</ProteinInfoList>");
        writer.flush();
    }


    public static void addAtoms(Protein model, CommandList commands) throws IOException{
        Command command = commands.firstWord(KeyWords.PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        BondParametersList bondParametersList  = new BondParametersList(parametersDirectory+"/" + MeshiPotential.BOND_PARAMETERS);
        AngleParametersList angleParametersList = new AngleParametersList(parametersDirectory+"/" + MeshiPotential.ANGLE_PARAMETERS);
        PlaneParametersList planeParametersList = new PlaneParametersList(parametersDirectory+"/" + MeshiPotential.PLANE_PARAMETERS);

        boolean hydrogenFailure = false;
        for (Atom atom : model.atoms()) {
            if (atom.isHydrogen() && atom.nowhere()) {
                PutHposLog puthLog = PutHpos.pos(atom, bondParametersList, angleParametersList);
                if (puthLog == null ) hydrogenFailure = true;
            }
        }
        AtomFinding.finalFindAtoms(model,bondParametersList,angleParametersList,planeParametersList);

    }
    private static ResidueList getConservedResidues(Protein model, CommandList commands) {
        ResidueList out = new ResidueList();
        String line = commands.firstWord("conserved").secondWord();
        String[] words = line.split(",");
        for (String word : words) {
            int number = Integer.valueOf(word);
            for (Residue residue : model.residues()) {
                if (residue.number() == number) {
                    out.add(residue);
                    Utils.println(residue + " marked as conserved");
                }
            }
        }
        return out;
    }

    private static void refreshBeta(RamachandranEnergy ramachandranEnergy) {
        ResidueTorsions rt1 = null, rt2 = null;
        ResidueTorsionsList residueTorsionsList = ramachandranEnergy.residueTorsionsList();
        for (ResidueTorsions residueTorsions : residueTorsionsList) {
            rt2 = residueTorsions;
            if (rt1 != null) {
                if (rt1.isBeta() && rt2.isBeta()) {
                    Residue residue = rt2.residue();
                    SecondaryStructure secondaryStructure = residue.getSecondaryStructure();
                    if (!secondaryStructure.equals(SecondaryStructure.SHEET))
                        Utils.println("Changing the secondary structure of " + residue + " to SHEET.");
                    rt2.residue().setSecondaryStructure(SecondaryStructure.SHEET);
                }
            }
            rt1 = rt2;
        }
    }


    private static class EvaluateLogger extends ProteinInfoListOld {
        String outFileName, infoFileName;
        MeshiWriter output;
        MeshiInfoXMLwriter info;
        ModelAnalyzer analyzer;
        String parentString;
        TotalEnergy energy;

        public EvaluateLogger(Protein model, String outFileName, TotalEnergy energy) {
            super("Optimization history of " + model);
            Protein nativeStructure;
            analyzer = new ModelAnalyzer(model, null, originalModel, energy, ResidueAlignmentMethod.IDENTITY);
            this.outFileName = outFileName;
            this.parentString = parentString;
            infoFileName = outFileName + ".xml";
            this.energy = energy;
        }

        public void setEnergy(TotalEnergy energy) {
            this.energy = energy;
            analyzer.setEnergy(energy);
        }


        public void log(ModelAnalyzer analyzer, String comment, boolean printFlag) throws AlignmentException{

            add(analyzer.analyze(comment));
            if (printFlag) {
                try {
                    output = new MeshiWriter(outFileName);
                    info = new MeshiInfoXMLwriter(infoFileName);
                    output.println(parentString);
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
                try {
                    print(output);
                } catch (IOException ex) {
                    Utils.throwException(this,ex,"failed to print");
                }
                try {
                    print(info);
                } catch (IOException ex) {
                    Utils.throwException(this,ex,"failed to print");
                }
                energy.on();
                try {
                    energy.evaluateAtoms();
                } catch (UpdateableException ex) {
                    System.out.println("log failed due to " + ex);
                    ex.printStackTrace();
                    throw new RuntimeException("quiting");
                } catch (EvaluationException ee) {
                    System.out.println("log failed due to " + ee);
                    ee.printStackTrace();
                    throw new RuntimeException("quiting");
                }
                Utils.colorByEnergy(model.atoms());
                analyzer.model.atoms().print(output);
                output.close();
            }
        }
         public void logXML(String comment,
                            boolean printFlag) throws IOException,EvaluationException,AlignmentException {
                add(analyzer.analyze(comment));
            if (printFlag) {
                try {
                    output = new MeshiWriter(outFileName);
                    info = new MeshiInfoXMLwriter(infoFileName);
                    output.println(parentString);
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
                print(output);
                print(info);
                energy.on();
                try {
                    energy.evaluateAtoms();
                } catch (UpdateableException ex) {
                    System.out.println("log failed due to " + ex);
                    ex.printStackTrace();
                    throw new RuntimeException("quiting");
                }
                Utils.colorByEnergy(model.atoms());
                analyzer.model.atoms().print(output);
                output.close();
            }
        }


        public void mcm(ArrayList<Score> scoreFunctions, TotalEnergy energy, String label) throws AlignmentException {
            analyzer.setEnergy(energy);
             long time = (new Date()).getTime() - energy.startTime;
            MeshiInfo infoList = energy.energyInfo();
            for (Score scoreFunction : scoreFunctions) {
                Utils.println(" Calculating score "+scoreFunction);
//                ((CombinedEnergyScore) scoreFunction).setEnergy(energy);
//                MeshiInfo lengthElement = new MeshiInfo(InfoType.LENGTH, InfoType.LENGTH.tag);
//                lengthElement.setValue(new Double(model.chain().numberOfNonDummyResidues()));
//                infoList.getChildren().add(lengthElement);

                double rmsFromOriginal;
                try {
                    rmsFromOriginal = Rms.rms(originalModel, model,ResidueAlignmentMethod.IDENTITY);
                }
                catch (Exception ex) {
                    ex.printStackTrace();
                    throw new RuntimeException(ex);
                }
                ((CombinedEnergyScore) scoreFunction).setChangeElement(rmsFromOriginal);
                MeshiInfo scores = scoreFunction.score(infoList);
                for (MeshiInfo score : scores.flatten())
                    label = label + "\" "+scoreFunction.toString()+"_"+score.type.tag+"=\""+score.getValue()+" ";
            }
        }

        public void mcm(ArrayList<Score> scoreFunctions,
                        TotalEnergy energy,
                        int i,
                        MCM.mcmStepResult mcmStepResult)  throws AlignmentException{
 //           mcm(optimizationScore, selectionScore, energy, "MCM\"  step=\""+i+"\" score=\""+score.score()+"\" result=\""+mcmStepResult+"\"  lastSuccess=\""+mcmStepResult.lastSuccess());
            mcm(scoreFunctions, energy, "MCM\"  step=\""+i+"\" result=\""+
                mcmStepResult+"\"  lastSuccess=\""+mcmStepResult.lastSuccess());
         }
    }


    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName", "startGdt_ha" ,"seed"};
        String[] arguments = getArguments(keys, argv);

        seed = Integer.parseInt(arguments[3]);
        System.out.println("seed " + seed);
        initRandom(seed);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        dsspFile   = inFileName+".dssp";
        datFile   = inFileName.substring(0,inFileName.length()-3)+"dat";
        startGdt_ha = Double.parseDouble(arguments[2]);
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
    }




}
