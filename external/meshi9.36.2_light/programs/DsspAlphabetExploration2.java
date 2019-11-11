/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.conservationContacts.ConservationContactsCreator;
import meshi.energy.conservationContacts.ConservationContactsHrCreator;
import meshi.energy.contacts.ContactsCreator;
import meshi.energy.goap.GoapCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZ5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZStd5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rapdf.RapdfCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.secondaryStructureAlphabet.FullSsAlphabetPredictionFeaturesCreator;
import meshi.energy.secondaryStructureAlphabet.SecondaryStructureAlphabetFeaturesCreator;
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
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore.RamachandranCoreCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZStdRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateType;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneEnergy;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.solvation.AtomEnvironmentCreator;
import meshi.energy.solvation.AtomEnvironmentCreatorSS;
import meshi.energy.solvation.SolvationCreator;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.geometry.ArcCos;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.*;
import meshi.parameters.MeshiPotential;
import meshi.parameters.SecondaryStructure;
import meshi.scoringFunctions.CombinedEnergyScore;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;
import meshi.util.info.MeshiInfoXMLwriter;
import meshi.util.info.ProteinInfoListOld;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

//import meshi.energy.one.OneCreator;

public class DsspAlphabetExploration2 extends MeshiProgram implements KeyWords {
    private enum TETHER_FLAG {
        ON, OFF, RESET_ON, RESET_OFF
    }

    private enum SUMMA_COOP_FLAG {
        ON, OFF
    }



    private enum RAMACH_COOP_FLAG {
        ON, OFF
    }



    private enum PROP_COOP_FLAG {
        ON, OFF
    }



    private enum RG_FLAG {
        ON, OFF
    }

    private enum SOLVATE_FLAG{
         ON, OFF
    }


    public static final String NAME = "Optimize";

    private static MCM mcm;
    private static Relaxation relaxation;
    private static LBFGS lbfgs;
    private static String inFileName, nativeFileName, outFileName,dsspFile;
    private static Protein model, originalModel;
    private static Boolean OK;
    private static CommandList commands;
    private static TotalEnergy minimizationEnergy;
    private static ArrayList<Score> scoreFunctions;
    private static TotalEnergy perturbationEnergy1, perturbationEnergy2, perturbationEnergy3, perturbationEnergy4;
    private static DistanceMatrix distanceMatrix;
    private static OptimizeLogger log;
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
    private static InflateCreator inflateByOtherModelCreator = new InflateCreator(InflateType.BY_OTHER_MODEL, new File("."), new MyFileFilter());
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
            //new DistanceConstraintsCreator(),
            new FullSsAlphabetPredictionFeaturesCreator(),
            new SecondaryStructureAlphabetFeaturesCreator()
    };
    private static EnergyCreator[] excludeFromMinimization = {};

    private static String parentString;

    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException, AlignmentException {
        init(argv);
        ArcCos.useFastArcCos();
        boolean stepwiseFlag = commands.keyExists("stepwise");
        // Chemistry
        model = null;

            parentString = getParentString(inFileName);
            model = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
            Utils.alignToX(model);
            originalModel = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
            addAtoms(model, commands);
            model.resetBonds();

            //model.getAtoms().backbone().freeze("Freezing backbon to allow side-chain relaxation");
            EnergyCreator[] relax1Creators = {new BondCreator(), new AngleCreator(), tetherAllCreator, residueTetherCreator};
            TotalEnergy firstRelaxEnergy = new TotalEnergy(model, relax1Creators, commands,"relaxEnergy");
            SteepestDecent steepestDecent = new SteepestDecent(firstRelaxEnergy,1,1000,200,0.00000001,0.5,1.5);
            try {
                steepestDecent.run();
            } catch (Exception ex) {
                throw new RuntimeException(ex);
            }
            Utils.println("nAtoms " + model.atoms().size());
            Utils.println("nResidues " + model.residues().size());
            //if (!OK) throw new RuntimeException("Too many atoms are missing. Cannot refine");
            ((MyFileFilter) inflateByOtherModelCreator.filter).reset(model);
            model.atoms().defrost();
            for (Atom atom : model.atoms()) {
                if (!atom.nowhere()) {
                    for (Atom neighbor : atom.bonded()) {
                        if (neighbor.nowhere()) atom.freeze();
                    }
                }
            }
    //      Secondary structure;
            if (!dsspFile.equals("NONE"))
                Utils.AssignFullDSSP(model, dsspFile);
            Utils.setSS(model, commands);
            Utils.println("secondary structure 1: ");
            int i = 0;
            for (Residue residue : model.residues()) {
                Utils.print(residue.getSecondaryStructure() + "  ");
                System.out.print(residue.type().nameOneLetter());//getSecondaryStructure().getOneLetterName());
                if (i++ % 10 == 0) Utils.println();
            }
        System.out.println();
            Utils.println();

            // conserved residues
            ResidueList conservedResidues = getConservedResidues(model, commands);

            // Work power

            minimizationEnergy = new TotalEnergy(model, excludeEnergy(energyCreators, excludeFromMinimization), commands,"minimization energy");
        //Instantiation TotalEnergy object may kill the other ones. Here we want all of them to be simultaneously alive (though not used).

        minimizationEnergy.terminator.reset();
            ((TetherEnergy) residueTetherCreator.term()).setResidue(commands);
            scoreFunctions = null;
            ArrayList<Score> energyScores = new ArrayList<Score>();



            log = new OptimizeLogger(model, originalModel, nativeFileName, outFileName, minimizationEnergy, parentString);
            log.mcm(scoreFunctions, minimizationEnergy, "BEGINNING");
            log.mcm(scoreFunctions, minimizationEnergy, "MCM_END\" step=\"" + 0);
            //log.log("BEGINNING", true);

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

    private static ArrayList<RamachandranEnergyElement> ramachEnergyPerResidue(RamachandranEnergy ramach) {
        ArrayList<RamachandranEnergyElement> out = new ArrayList<RamachandranEnergyElement>();
        for (Object o : ramach.elementsList()) {
            RamachandranEnergyElement element = (RamachandranEnergyElement) o;
            double e = element.evaluate() / element.weight();
            if (e > 5) Utils.println("Ramachandran element " + element + " ; minimizationEnergy = " + e);
            if (e > RamachandranEnergyElement.THRESHOLD) out.add(element);
        }
        return out;
    }


    private static MCM getMCM(Protein model,
                                  TotalEnergy minimizationEnergy,
                                  ArrayList<Score> scoreFunctions,
                                  TotalEnergy perturbationEnergy1,
                                  TotalEnergy perturbationEnergy2,
                                  TotalEnergy perturbationEnergy3,
                                  TotalEnergy perturbationEnergy4,
                                  CommandList commands,
                                  ResidueList conservedResidues) throws UpdateableException,EvaluationException {
        Score optimizationScore = null;
        for (Score score:scoreFunctions)
            if (score.toString().equals("optimizationScore"))  {
                optimizationScore = score;
                break;
            }
        if (optimizationScore == null)
            throw new RuntimeException("No optimization score.");

        return getMCM(model,minimizationEnergy,scoreFunctions,optimizationScore,
                                         perturbationEnergy1,perturbationEnergy2,perturbationEnergy3,perturbationEnergy4,
                                         commands, conservedResidues, MCM.mcmMode.OPTIMIZATION);
        }
 /*   private static Relaxation getRelaxation(Protein model,
                                  TotalEnergy minimizationEnergy,
                                  ArrayList<Score> scoreFunctions,
                                  TotalEnergy perturbationEnergy1,
                                  TotalEnergy perturbationEnergy2,
                                  TotalEnergy perturbationEnergy3,
                                  TotalEnergy perturbationEnergy4,
                                  CommandList commands,
                                  ResidueList conservedResidues) throws UpdateableException,EvaluationException {
            return (Relaxation) getMCM(model, minimizationEnergy, scoreFunctions,
                                                                  perturbationEnergy1, perturbationEnergy2,perturbationEnergy3,perturbationEnergy4,
                                                                  commands, conservedResidues, MCM.mcmMode.RELAXATION);
        }
   */
    private static MCM getMCM(Protein model,
                              TotalEnergy minimizationEnergy,
                              ArrayList<Score> scoreFunctions,
                              Score optimizationScore,
                              TotalEnergy perturbationEnergy1,
                              TotalEnergy perturbationEnergy2,
                              TotalEnergy perturbationEnergy3,
                              TotalEnergy perturbationEnergy4,
                              CommandList commands,
                              ResidueList conservedResidues, MCM.mcmMode mode) throws UpdateableException,EvaluationException {
        Perturbation perturbation1 = new PerturbationByMinimization(perturbationEnergy1, commands, model, conservedResidues, " perturbation with a complete energy function");
        Perturbation perturbation2 = new PerturbationByMinimization(perturbationEnergy2, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation perturbation3 = new PerturbationByMinimization(perturbationEnergy3, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation perturbation4 = new PerturbationByMinimization(perturbationEnergy4, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation[] perturbations = {perturbation1, perturbation2, perturbation3};
        Perturbation perturbation = new CombinedPerturbation(perturbations);
        //new ScmodPerturbations(model, commands,conservedResidues)};
        minimizationEnergy.setCurrent();
        Minimizer minimizer;
        minimizer = Utils.getLBFGS(minimizationEnergy, commands, MINIMIZE);
        double initialTemperature = commands.firstWordFilter(MC_MINIMIZATION).secondWord(INITIAL_TEMPERATURE).thirdWordDouble();
        double finalTemperature = commands.firstWordFilter(MC_MINIMIZATION).secondWord(FINAL_TEMPERATURE).thirdWordDouble();
        int nSteps = commands.firstWordFilter(MC_MINIMIZATION).secondWord(MAX_STEPS).thirdWordInt();
        TemperatureGenerator temperatureGenerator = new TemperatureGenerator(initialTemperature, finalTemperature, nSteps);
        //AbstractEnergy[] excludedTerms = {inflateCreator.term(), tetherAllCreator.term()};
        if (mode == MCM.mcmMode.RELAXATION )
            return new Relaxation(minimizationEnergy, scoreFunctions, optimizationScore, minimizer, perturbation, temperatureGenerator, nSteps);
        return new MCM(minimizationEnergy, scoreFunctions, optimizationScore, minimizer, perturbation, temperatureGenerator, nSteps);
    }

    private static Optimizer.OptimizerStatus minimize(Minimizer lbfgs, TotalEnergy energy,
                                 TETHER_FLAG tetherFlag,
                                 RG_FLAG rgFlag,SOLVATE_FLAG solvateFlag) throws OptimizerException {

        return minimize(lbfgs, energy, RAMACH_COOP_FLAG.ON,
                PROP_COOP_FLAG.ON,
                SUMMA_COOP_FLAG.ON,
                tetherFlag,
                rgFlag,solvateFlag);
    }

    private static Optimizer.OptimizerStatus minimize(Minimizer lbfgs,TotalEnergy energy,
                                 RAMACH_COOP_FLAG rcFlag,
                                 PROP_COOP_FLAG prFlag,
                                 SUMMA_COOP_FLAG scFlag,
                                 TETHER_FLAG tetherFlag,
                                 RG_FLAG rgFlag,
                                 SOLVATE_FLAG solvateFlag)  throws OptimizerException{


        if ((rcFlag == RAMACH_COOP_FLAG.OFF) &&
            (cooperativeZRamachandranCreator.term() != null)) {
            cooperativeZRamachandranCreator.term().off();
            cooperativeZStdRamachandranCreator.term().off();
        }
        if ((prFlag == PROP_COOP_FLAG.OFF) &&
            (cooperativeZPropensityCreator.term() != null)){
            cooperativeZPropensityCreator.term().off();
            cooperativeZStdPropensityCreator.term().off();
        }
        if ((scFlag == SUMMA_COOP_FLAG.OFF) &&
            (cooperativeZSummaCreator.term()!= null)) {
            cooperativeZSummaCreator.term().off();
        }
        if (tetherFlag == TETHER_FLAG.OFF) {
            tetherAllCreator.term().off();
        } else if (tetherFlag == TETHER_FLAG.RESET_OFF) {
            ((TetherEnergy) tetherAllCreator.term()).reset();
            tetherAllCreator.term().off();
        } else if (tetherFlag == TETHER_FLAG.RESET_ON) {
            ((TetherEnergy) tetherAllCreator.term()).reset();
        }
        if (solvateFlag == SOLVATE_FLAG.OFF) {
            solvationCreator.term().off();
        }
        RamachandranEnergy.resetRamachandran(energy);

        try {
            return lbfgs.run();
        } catch (Exception ex) {
            writeFailureXml(model, ex);
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }



    private static class OptimizeLogger extends ProteinInfoListOld implements Logger {
        String outFileName, infoFileName;
        MeshiWriter output;
        MeshiInfoXMLwriter info;
        ModelAnalyzer analyzer;
        String parentString;
        TotalEnergy energy;

        public OptimizeLogger(Protein model, Protein originalModel, String nativeFileName, String outFileName, TotalEnergy energy, String parentString) {
            super("Optimization history of " + model);
            Protein nativeStructure;
            if (!(nativeFileName.equals("NONE") ||nativeFileName.equals("none"))) {
                nativeStructure = Utils.getProtein(commands, nativeFileName, ResidueExtendedAtoms.creator, Utils.defaultExceptionHandler);
            } else {
                nativeStructure = null;
            }
            analyzer = new ModelAnalyzer(model, nativeStructure, originalModel, energy, ResidueAlignmentMethod.IDENTITY);
            this.outFileName = outFileName;
            this.parentString = parentString;
            infoFileName = outFileName + ".xml";
            this.energy = energy;
        }

        public void setEnergy(TotalEnergy energy) {
            this.energy = energy;
            analyzer.setEnergy(energy);
        }

        public double rms() {
            if(nativeFileName.equals("NONE") ||nativeFileName.equals("none") ) return -1;
            else {
                try{
                    return analyzer.rms();
                }catch (Exception ex ) {
                    writeFailureXml(model, ex);
                    throw new RuntimeException(ex.getMessage());}
            }
        }

        public void log(String comment, boolean printFlag) throws AlignmentException{

            add(analyzer.analyze(comment));
            if (printFlag) {
                try {
                    output = new MeshiWriter(outFileName);
                    info = new MeshiInfoXMLwriter(infoFileName);
                    output.println(parentString);
                } catch (Exception ex) {
                    writeFailureXml(model, ex);
                    throw new RuntimeException(ex);
                }
                try {
                    print(output);
                } catch (IOException ex) {
                    writeFailureXml(model, ex);
                    Utils.throwException(this,ex,"failed to print");
                }
                try {
                    print(info);
                } catch (IOException ex) {
                    writeFailureXml(model, ex);
                    Utils.throwException(this,ex,"failed to print");
                }
                energy.on();
                try {
                    energy.evaluateAtoms();
                } catch (UpdateableException ex) {
                    writeFailureXml(model, ex);
                    System.out.println("log failed due to " + ex);
                    ex.printStackTrace();
                    throw new RuntimeException("quiting");
                } catch (EvaluationException ee) {
                    writeFailureXml(model, ee);
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
                    writeFailureXml(model, ex);
                    throw new RuntimeException(ex);
                }
                print(output);
                print(info);
                energy.on();
                try {
                    energy.evaluateAtoms();
                } catch (UpdateableException ex) {
                    writeFailureXml(model, ex);
                    System.out.println("log failed due to " + ex);
                    ex.printStackTrace();
                    throw new RuntimeException("quiting");
                }
                Utils.colorByEnergy(model.atoms());
                analyzer.model.atoms().print(output);
                output.close();
            }
        }


        public void mcm(ArrayList<Score> scoreFunctions, TotalEnergy energy, String label)  throws AlignmentException{
            analyzer.setEnergy(energy);
             long time = (new Date()).getTime() - energy.startTime;
            MeshiInfo infoList = energy.energyInfo();
            if (scoreFunctions != null) {
                for (Score scoreFunction : scoreFunctions) {
                    Utils.println(" Calculating score " + scoreFunction);
                    //                ((CombinedEnergyScore) scoreFunction).setEnergy(energy);
                    //                MeshiInfo lengthElement = new MeshiInfo(InfoType.LENGTH, InfoType.LENGTH.tag);
                    //                lengthElement.setValue(new Double(model.chain().numberOfNonDummyResidues()));
                    //                infoList.getChildren().add(lengthElement);

                    double rmsFromOriginal;
                    try {
                        rmsFromOriginal = Rms.rms(originalModel, model, ResidueAlignmentMethod.IDENTITY);
                    } catch (Exception ex) {
                        writeFailureXml(model, ex);
                        ex.printStackTrace();
                        throw new RuntimeException(ex);
                    }
                    ((CombinedEnergyScore) scoreFunction).setChangeElement(rmsFromOriginal);
                    MeshiInfo scores = scoreFunction.score(infoList);
                    for (MeshiInfo score : scores.flatten())
                        label = label + "\" " + scoreFunction.toString() + "_" + score.type.tag + "=\"" + score.getValue() + " ";
                }
            }//if
 	        log(label+"\" time=\""+time, true);
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
        String[] keys = {"commands", "inFileName", "dsspFile", "nativeStructure", "outFile", "seed"};
        String[] arguments = getArguments(keys, argv);

        seed = Integer.parseInt(arguments[5]);
        System.out.println("seed " + seed);
        initRandom(seed);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        dsspFile   = arguments[2];
        nativeFileName = arguments[3];
        outFileName = arguments[4];
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
    }

    private static String getParentString(String fileName) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String out;
        while ((out = reader.readLine()) != null) {
            if (out.startsWith("PARENT"))
                return out;
        }
        return "PARENT N/A";
    }


    private static EnergyCreator[] excludeEnergy(EnergyCreator[] source, EnergyCreator[] exclude) {
        EnergyCreator[] out = new EnergyCreator[source.length - exclude.length];
        int j = 0;
        for (int i = 0; i < source.length; i++) {
            if (notFound(source[i], exclude)) {
                out[j] = source[i];
                j++;
            }
        }
        return out;
    }

    private static boolean notFound(EnergyCreator energyCreator, EnergyCreator[] energyCreators) {
        for (int i = 0; i < energyCreators.length; i++) {
            if (energyCreators[i].equals(energyCreator)) return false;
        }
        return true;
    }

    private static class MyFileFilter implements Filter {
        private String prefix = null;
        Protein thisModel;

        public void reset(Protein thisModel) {
            this.thisModel = thisModel;
            int index = thisModel.name().indexOf('.');
            if (index != -1) prefix = thisModel.name().substring(0, index);
            else prefix = thisModel.name();

        }

        public boolean accept(Object obj) {
            if (prefix == null) throw new RuntimeException("prefix is " + prefix);
            File file = (File) obj;
            String path = file.getAbsolutePath();
            if (!path.endsWith("pdb")) return false;
            if (path.indexOf("out") == -1) return false;
            if (file.getName().startsWith(prefix)) {
                try {
                    double rms = Rms.rms(thisModel, Protein.getCAproteinFromApdbFile(file), ResidueAlignmentMethod.IDENTITY);
                    if (rms < 1) return false;
                } catch (Exception ex ) {
                    writeFailureXml(model, ex);
                    throw new RuntimeException(ex.getMessage());}
                return true;
            }
            return false;
        }
    }

    private static void writeFailureXml(Protein model, Exception exception)  {
        MeshiWriter writer;
        try {
            writer = new MeshiWriter(outFileName + ".xml");
        } catch (IOException ex) {throw  new RuntimeException("Cannot write failure XML file after exception:\n"+exception+"\n"+"Due to "+ex);}
        writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\" ?> ");
        writer.println("<ProteinInfoList name=\"Failure report for Protein: " + model.sourceFile() + "\">");
        writer.print("<ProteinInfo  value=\"MCM_END\" step=\"0\" ");
//        for (Score scoreFunction : scoreFunctions) {
//            writer.print(scoreFunction.toString()+"_weightedMedianScore=\"0.05\" "+scoreFunction.toString()+"_interdecile=\"0\" ");
//        }
        writer.println("time=\"0\" fileName=\""+inFileName+"\" >");
        writer.println("<Exception>\n"+ exception + "\n</Exception>" );
        writer.println("</ProteinInfo>\n" + "</ProteinInfoList>\n");
        writer.close();
    }
}
