/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 17/09/2009
 * Time: 14:12:03
 * To change this template use File | Settings | File Templates.
 */
package programs;

import meshi.PDB.PdbLineATOM;
import meshi.applications.helixParabolaAttractor.HelixParabolaAttractorCreator;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.beta.BetaCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsEnergy;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.rg.RgEnergy;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.fragments.DresserFullFrag;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.BackboneResidue;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.LBFGS;
import meshi.parameters.MeshiPotential;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.ResidueSequence;
import meshi.sequences.SequenceAlignmentCell;
import meshi.sequences.SequenceAlignmentColumn;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.string.StringList;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;


public class AbInitio extends MeshiProgram implements MeshiPotential, KeyWords {
    public static final double HIGHEST_RAMACH_ELEMENT_FACTOR = 100;
    private static CommandList commands;
    private static int length;
    private static String modelFileName;
    private static boolean debug = true;
    public static final String NAME = "AbInitio";
    private static Protein reference;

    private static final HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
    private static final AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
    private static final RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
    private static final FlatRamachCreator flatRamahCreator = new FlatRamachCreator();
    private static final RgCreator rgCreator = new RgCreator();
    private static final BetaCreator betaCreator = new BetaCreator();
    private static final HelixParabolaAttractorCreator helixCreator = new HelixParabolaAttractorCreator();
    private static final HydrogenBondsPairsCreator hydrogenBondPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);

    public static void main(String[] args) throws Exception {
        System.out.print("command line: ");
        for (String arg : args) System.out.print(arg + " ");
        System.out.println();
        int seed = Integer.parseInt(args[4]);
        commands = Utils.init(args, 5, seed, "Usage: java -XmxNNNm " + NAME + "  <commands file name><sequnce file name> <native pdb file name> <output file name> <seed>\n\n" +
                "NNN is the size of the expected memory requirement in MegaBytes.");
        MolecularSystem molecularSystem = new MolecularSystem();
        String inFileName, referenceFileName, outFileName, referenceOutFileName;
        Protein reference;
        AbInitioLogger log;
        ResidueSequence targetSequence;
        TotalEnergy energy;
        int numberOfMinimizations;
        DistanceMatrix distanceMatrix;

        EnergyCreator[] firstEnergyCreators = {
                new BondCreator(),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                new AtomicPairwisePMFSummaCreator(),
                new FlatRamachCreator(),
                new HelixParabolaAttractorCreator()
        };

        EnergyCreator[] energyCreators = {
                ramachCreator,
                summaCreator,
                rgCreator,
                betaCreator,
                hydrogenBondsCreator,
                hydrogenBondPairsCreator,
                new BondCreator(),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                helixCreator,
                flatRamahCreator,
                new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
                new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
        };


        inFileName = args[1];
        referenceFileName = args[2];
        outFileName = args[3];


        referenceOutFileName = referenceFileName + ".out";
        reference = Utils.getProtein(commands, referenceFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
        Utils.AssignDSSP(reference, commands, KeyWords.SECONDARY_STRUCTURE);
        distanceMatrix = reference.atoms().molecularSystem().getDistanceMatrix();
        energy = new TotalEnergy(reference, energyCreators, commands,"generic energy (abInitio 1");
        energy.evaluate();
        MeshiWriter writer = new MeshiWriter(referenceOutFileName);
        String referenceEnergy = energy.report(0);
        writer.println(referenceEnergy);
        writer.close();

        RgCalculator rgCalculator = new RgCalculator(reference.atoms());
        String name = reference.name();
        double radius2 = reference.atoms().radius2();

        numberOfMinimizations = commands.firstWord("numberOfMinimizations").secondWordInt();
        if (numberOfMinimizations < 3)
            throw new RuntimeException("numberOfMinimizations = " + numberOfMinimizations + " must be larger than 3.");

        System.out.println("reference radius : " + Math.sqrt(radius2) + " , " + rgCalculator.rg());

        Command command = commands.firstWord(PARAMETERS_DIRECTORY);
        String FragmentsFileName = command.secondWord();
        command = commands.firstWord(DRESSER_FRAGMENTS);
        FragmentsFileName = FragmentsFileName + "/" + command.secondWord();
        DresserFullFrag dresser = new DresserFullFrag(FragmentsFileName);

        targetSequence = getResidueSequence(inFileName, name);
        double[] bais = {1, 0.5, 0.5};
        Protein model;
        model = createInitialModel(targetSequence, bais, dresser, firstEnergyCreators,molecularSystem);
        if (model == null) throw new RuntimeException("Failed to create an initial model.");
        distanceMatrix = model.atoms().molecularSystem().getDistanceMatrix();
        energy = new TotalEnergy(model, energyCreators, commands,"generic energy (abInitio 2");
        log = new AbInitioLogger(model, reference, energy, commands, outFileName);
        model = refine(model, energy, log, numberOfMinimizations);
        if (model == null) throw new RuntimeException("Failed to refine the initial model.");
        log.printModel(numberOfMinimizations - 1);

    }//main

    public static Protein createInitialModel(ResidueSequence targetSequence, double[] bais,
                                             DresserFullFrag dresser,
                                             EnergyCreator[] energyCreators,MolecularSystem molecularSystem) throws UpdateableException, EvaluationException{

        System.out.println("use  assignRandomCaCoordinatesArounASphere");
        new MolecularSystem();
        System.out.println("Set the begenning position !");
        Protein backboneModel = new Protein(targetSequence, name, new BackboneResidue("creator"));
        int n = backboneModel.atoms().size();
        double targetRG = 1.8 * Math.exp(0.4 * Math.log(2.0 * n) - 0.25);
        System.out.println("Target RG for random chain " + targetRG);
        double rg = 999;
        while (rg > targetRG) {
            Utils.assignRandomCaCoordinatesInACubeRandomLength(backboneModel.chain(),
                    backboneModel.chain().numberOfNonDummyResidues() / 8,
                    backboneModel.chain().numberOfNonDummyResidues() / 4, bais);
            rg = Utils.radiusOfGyration(backboneModel.atoms());
            System.out.println("RG = " + rg + "  targetRG = " + targetRG);
        }
        AtomList tempAtomList = dresser.dressProt(backboneModel.atoms().CAFilter());
        Utils.assignBackboneCoordinates(backboneModel.atoms(), tempAtomList);
        Protein model = new Protein(backboneModel.name());
        ResidueExtendedAtomsCreator creator = new ResidueExtendedAtomsCreator();
        new MolecularSystem();
        Chain backboneChain = backboneModel.chains().get(0);
        Chain chain = new Chain(backboneChain.name(), Chain.GENERIC_CHAIN_NUMBER, model);
        for (Object aBackboneChain : backboneChain) {
            Residue residue = (Residue) aBackboneChain;
            if (residue.dummy()) chain.add(residue);
            else {
                Residue newResidue = creator.create(residue.getAtoms(), residue.ID(), residue.mode,molecularSystem);
                newResidue.setSecondaryStructure(residue.getSecondaryStructure());
                chain.add(newResidue);
            }
        }
        model.addChain(chain);
        Scmod.scmod(commands, model, 1, 50);
        DistanceMatrix distanceMatrix = model.atoms().molecularSystem().getDistanceMatrix();
        TotalEnergy firstEnergy = new TotalEnergy(model, energyCreators, commands,"generic energy (createInitialModel)");
        try {
            LBFGS lbfgs = new LBFGS(firstEnergy, 0.02, 2000, 500);
            lbfgs.run();
        } catch (Exception ex) {
            return null;
        }
        return model;
    }

    private static Protein refine(Protein model, TotalEnergy energy, AbInitioLogger log,
                                  int numberOfMinimizations) throws UpdateableException, EvaluationException {
        boolean headerFlag = true;
        for (int i = 0; i < numberOfMinimizations; i++) {
            System.out.println("Minimization # " + i + " of " + numberOfMinimizations);
            energy.setCurrent();
            LBFGS lbfgs = Utils.getLBFGS(energy, commands, RELAX);
            if (i <= 1) {            //First minimization. collapse with RG . The beta term is not used to prevent bias of the secondary structure towards local interactions
                betaCreator.term().off();
                rgCreator.term().off();
                ((HydrogenBondsPairsEnergy) hydrogenBondPairsCreator.term()).off();
                ((RamachandranSidechainEnergy) ramachCreator.term()).scaleWeight(10);
            } else if (i == 2) {
                betaCreator.term().off();
                rgCreator.term().on();
                ((RgEnergy) rgCreator.term()).scaleWeight(0.1);
                ((RamachandranSidechainEnergy) ramachCreator.term()).scaleWeight(0.1);
                ((HydrogenBondsPairsEnergy) hydrogenBondPairsCreator.term()).on();
            } else if (i == 3) {
                betaCreator.term().off();
                rgCreator.term().off();
                ((RamachandranSidechainEnergy) ramachCreator.term()).scaleWeight(0.1);
            } else if (i == 4) {
                betaCreator.term().off();
                rgCreator.term().on();
                ((RgEnergy) rgCreator.term()).scaleWeight(5);        // to get weight of 0.5   of the assigned by the user
            } else if (i == 5) {
                betaCreator.term().off();
                rgCreator.term().off();
                ((RgEnergy) rgCreator.term()).scaleWeight(2);        // to get weight as assigned by the user
            } else if (i > 5) {
                helixCreator.term().off();
                flatRamahCreator.term().off();
                if (i % 2 == 1) {
                    betaCreator.term().off();
                    rgCreator.term().off();
                } else {
                    betaCreator.term().on();
                    //((BetaEnergy)betaCreator.term()).scale(1.2);
                    rgCreator.term().on();
                }
            }
            try {
                lbfgs.run();
                log.addLine(i, "Initialization Try0 ", headerFlag);
                headerFlag = false;
                if (i < (numberOfMinimizations - 2))
                    Scmod.scmod(commands, model, 4, 50);
            }
            catch (Exception ex) {
                System.out.println("Cannot continue. Caught \n" + ex);
                return null;
            }
        }
        return model;
    }


    private static ResidueSequence getResidueSequence(String fileName, String name) {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        StringList lines = new StringList(reader);
        ResidueSequence out = new ResidueSequence(name);
        for (Object line1 : lines) {
            String line = (String) line1;
            StringTokenizer words = new StringTokenizer(line);
            String word1 = words.nextToken();
            String word2 = words.nextToken();
            String word3 = words.nextToken();
            int number = Integer.parseInt(word1);
            // char seqC = word2.charAt(0);
            char seqC = ResidueNameConvertor.three2one(word2);
            char ssC = word3.charAt(0);

            if ("ACDEFGHIKLMNPQRSTVWYZ".indexOf(seqC) < 0)
                throw new RuntimeException("weird seq char: " + seqC);
            if ("HE-".indexOf(ssC) < 0)
                throw new RuntimeException("weird ss  char: " + ssC);
            if (ssC == '-')
                ssC = 'C';

            SequenceAlignmentCell cell = new SequenceAlignmentCell(seqC, number);
            SequenceAlignmentColumn column = new SequenceAlignmentColumn(cell);
            SecondaryStructure ss = SecondaryStructure.secondaryStructure(ssC);
            cell.addAttribute(ss);
            out.add(column);
        }
        return out;
    }


    private static void summary(TotalEnergy e, int stage) {
        ArrayList energyTerms = e.energyTerms();
        ArrayList energyValues = e.energyValues();
        Iterator values = energyValues.iterator();

        Double value;
        System.out.println("\n\nSUMMARY " + stage + ":\n\tTotalEnergy: \t" + e.energy());
        for (Object energyTerm : energyTerms) {
            value = (Double) values.next();
            if ((energyTerm != null) && (value != null))
                System.out.println("\t" + ((AbstractEnergy) energyTerm).comment() + ": \t" + value);
        }
        System.out.println();
    }


    protected static CommandList init(String[] argv) throws Exception {
        //int nAtoms;
        //double initialRadius;
        CommandList commands;
        String help = "Usage: java HydrogenTest <commands file name> <modelFileName> <modelDsspFileName> seed";

        if ((argv.length == 0) || getFlag("-help", argv)) {
            throw new Exception(help);
        }

        try {
            //commands = new CommandList(getOrderedArgument(argv));
            String commandsName = getOrderedArgument(argv);
            String[] tmpArray = {commandsName};
            commands = Utils.init(tmpArray, "HydrogenTest");
            modelFileName = getOrderedArgument(argv);
            if (modelFileName == null) throw new RuntimeException(help);
            System.out.println("# initial model file name is " + modelFileName);
            reference = new Protein(modelFileName + ".pdb", new PdbLineATOM(), new ResidueExtendedAtomsCreator(), null);
            length = reference.chain().size() - reference.chain().firstNonDummyResidue().number();
            System.out.println("# length: " + length);
            initRandom(argv);
        }
        catch (Exception e) {
            throw new Exception(help + "\n" + e.toString());
        }
        return commands;
    }


    private static class AbInitioLogger extends StringList {
        private boolean referenceExists;
        private Protein model;
        private Protein reference;
        private String outFileName;
        private CommandList commands;
        private MeshiWriter output;
        TotalEnergy energy;

        public AbInitioLogger(Protein model, Protein reference, TotalEnergy energy, CommandList commands, String outFileName) {
            this.reference = reference;
            this.commands = commands;
            this.outFileName = outFileName;
            this.model = model;
            this.energy = energy;
        }

        public void addLine(int i, String label, boolean headerFlag) {
            /*   Utils.RmsGdtEnergy(this, model, reference, energy, label+i+"  ", headerFlag); */
        }

        public void printModel(int step) throws UpdateableException,EvaluationException {
            try {
                output = new MeshiWriter(outFileName);
            }
            catch (Exception ex) {
                throw new RuntimeException(ex);
            }
            print(output);
            energy.on();
            helixCreator.term().off();
            flatRamahCreator.term().off();
            if (step % 2 == 1) {
                betaCreator.term().off();
                rgCreator.term().off();
            }
            energy.evaluateAtoms();
            Utils.colorByEnergy(model.atoms());
            model.atoms().print(output);
            output.close();
        }
    }

    public static class EnergyAttribute implements MeshiAttribute {
        public final Double energy;

        public EnergyAttribute(Double energy) {
            this.energy = energy;
        }

        public int key() {
            return MeshiAttribute.ENERGY_ATTRIBUTE;
        }
    }

}
