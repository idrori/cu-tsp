package programs;

import meshi.energy.beta.*;
import meshi.energy.hydrogenBond.*;
import meshi.energy.hydrogenBondsAngle.*;
import meshi.energy.hydrogenBondsPairs.*;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.*;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.*;
import meshi.energy.rg.*;
import meshi.energy.simpleEnergyTerms.angle.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.*;
import meshi.energy.simpleEnergyTerms.outOfPlane.*;
import meshi.energy.simpleEnergyTerms.plane.*;
import meshi.energy.simpleEnergyTerms.tether.*;
import meshi.energy.solvation.*;
import meshi.molecularElements.Protein;
import meshi.molecularElements.extendedAtoms.*;
import meshi.optimizers.*;
import meshi.parameters.MeshiPotential;
import meshi.sequences.AlignmentException;
import meshi.util.*;
import meshi.energy.*;
import meshi.util.info.InfoType;

import java.io.*;

/**
 *
 */
public class ColorByEnergy extends MeshiProgram implements KeyWords {
    private static TetherCreator tetherCreator = new TetherCreator(InfoType.TETHER_ALL);
    private static HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
    private static SolvationCreator solvationCreator = new SolvationCreator(1.0, 1.0, 1.0, 1.0, 1.0);
    private static RgCreator rgCreator = new RgCreator();

    private static CompositePropensityCreator propensityCreator = new CompositePropensityCreator();
    private static CooperativeZPropensityCreator cooperativeZPropensityCreator = new CooperativeZPropensityCreator(propensityCreator);
    private static CooperativeZStdPropensityCreator cooperativeZStdPropensityCreator = new CooperativeZStdPropensityCreator(cooperativeZPropensityCreator);

    private static AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
    private static CooperativeZ5typesSummaCreator cooperativeZSummaCreator = new CooperativeZ5typesSummaCreator(summaCreator);

    private static RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
    private static CooperativeZRamachandranCreator cooperativeZRamachandranCreator = new CooperativeZRamachandranCreator(ramachCreator);
    private static CooperativeZStdRamachandranCreator cooperativeZStdRamachandranCreator = new CooperativeZStdRamachandranCreator(cooperativeZRamachandranCreator);
    private static RamachandranCreator ramachandranCreator = new RamachandranCreator();
    private static BetaCreator betaCreator = new BetaCreator();
    private static EnergyCreator[] energyCreators = {new BondCreator(),
            new AngleCreator(),
            new PlaneCreator(),
            new OutOfPlaneCreator(),
            ramachandranCreator,
            ramachCreator,
            cooperativeZRamachandranCreator,
            cooperativeZStdRamachandranCreator,
            propensityCreator,
            cooperativeZPropensityCreator,
            cooperativeZStdPropensityCreator,
            summaCreator,
            cooperativeZSummaCreator,
            solvationCreator,
            hydrogenBondsCreator,
            new HydrogenBondsPairsCreator(hydrogenBondsCreator),
            new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
            new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
            rgCreator,
            betaCreator,
            tetherCreator
    };
    private static int seed;
    public static final String NAME = "ColorByEnergy";
    private static String inFileName, nativeFileName, outFileName;

    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException,AlignmentException{
        CommandList commands = init(argv);

        Protein model = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
        boolean OK = Utils.addAtoms(model, false, commands,new PlaneParametersList(EnergyCreator.parametersDirectory(commands) +
                "/" + MeshiPotential.PLANE_PARAMETERS),false);
    }

    private static CommandList init(String[] argv) {
        if (argv.length >= 5) seed = Integer.parseInt(argv[4]);
        else seed = -1;
        System.out.println("seed " + seed);
        CommandList commands = Utils.init(argv, 5, seed, "Usage: java -XmxNNNm " + NAME + "  <commands file name><model pdb file name> <native pdb file name> <output file name> <seed> \n\n" +
                "NNN is the size of the expected memory requirement in MegaBytes.");
        inFileName = argv[1];
        nativeFileName = argv[2];
        outFileName = argv[3];
        return commands;
    }
}

