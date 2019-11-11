/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public interface KeyWords {
    //---------------------------- protein -------------------------------
    final Key SEQUENCE = new Key("sequence");
    final Key SECONDARY_STRUCTURE = new Key("secondary");
    final Key MODEL_NUMBER = new Key("modelNumber");
    final Key CLASH_DISTANCE = new Key("clashDistance");
    final Key MAX_CLASHES = new Key("maxNumberOfClashes");
    final Key N_TRIES = new Key("nTries");

    //---------------------------- optimization -------------------------------
    final Key OPTIMIZER = new Key("optimizer");
    final Key MINIMIZE = new Key("minimize"); // First word for full minimization
    final Key RELAX = new Key("relax"); // First word for short relaxation minimization with SteepestDecent
    final Key TOLERANCE = new Key("tolerance");
    final Key MAX_STEPS = new Key("maxSteps");
    final Key REPORT_EVERY = new Key("reportEvery");
    final Key RESTART_EVERY = new Key("restartEvery"); // for Conjugate Gradients
    final Key STEEPEST_DECENT = new Key("steepestDecent");
    final Key BFGS = new Key("bfgs");
    final Key LBFGS = new Key("lbfgs");
    final Key CG = new Key("cg");

    //--------------------------------- Function  prediction ----------------------
    final Key STRUCTURE_NAMES = new Key("structureNames");
    final Key CSAonly_FILES_LOCATION_PATH = new Key("csaOnlyFilesLocationPath");
    final Key DISTANCE_FROM_CENTROID_ENERGY = new Key("distanceFromCentroid");
    //--------------------------------- energy  -----------------------------
    final Key PARAMETERS_DIRECTORY = new Key("parameters");

    //--------------------------------- energy terms -----------------------------
    final Key LENNARD_JONES_CA = new Key("LennardJonesCa");
    final Key LENNARD_JONES = new Key("unsupportedEnergies");
    final Key ANGLE_ENERGY = new Key("angleEnergy");
    final Key BOND_ENERGY = new Key("bondEnergy");
    final Key CONSENSUS_ENERGY = new Key("consensusEnergy");
    final Key PLANE_ENERGY = new Key("planeEnergy");
    final Key OUT_OFPLANE_ENERGY = new Key("outOfPlaneEnergy");
    final Key TEMPLATE_DISTANCE_CONSTRAINTS = new Key("templateDistanceConstraints");
    final Key DISTANCE_CONSTRAINTS_ENERGY = new Key("distanceConstraintsEnergy");
    final Key INFLATE_ENERGY = new Key("inflateEnergy");
    final Key VOLUME_CONSTRAINT = new Key("volumeConstraintWeight");
    final Key HYDROGEN_BONDS = new Key("hydrogenBonds");
    final Key HYDROGEN_BONDS_PAIRS = new Key("hydrogenBondsPairs");
    final Key MAX_DISTANCE = new Key("maxDistance");
    final Key MAX_ANGLE = new Key("maxAngle");
    final Key HELIX_ATTRACTOR = new Key("helixAttractor");
    final Key BETA_ATTRACTOR = new Key("betaAttractor");
    final Key BETA = new Key("beta");


    final Key TWO_TORSIONS_ENERGY = new Key("twoTorsionsEnergy");
    final Key FLAT_RAMACH_ENERGY = new Key("flatRamachEnergy");
    final Key PROPENSITY_TORSION_ENERGY = new Key("propensityTorsionEnergy");
    final Key ALPHA_ANGLE_ENERGY = new Key("alphaAngleEnergy");
    final Key ALPHA_TORSION_ENERGY = new Key("alphaTorsionEnergy");
    final Key SOLVATE_ENERGY = new Key("solvateEnergy");
    final Key SOLVATE_ENERGY_SC_POLAR = new Key("solvateEnergySCpolar");
    final Key SOLVATE_ENERGY_SC_CARBON = new Key("solvateEnergySCcarbon");
    final Key SOLVATE_ENERGY_BB_POLAR = new Key("solvateEnergyBBpolar");
    final Key SOLVATE_ENERGY_HB = new Key("solvateEnergyHB");
    final Key EXCLUDED_VOL = new Key("excludedVolume");
    final Key ELECTROSTATICS = new Key("electrostatics");
    final Key DIELECTRIC_CONSTANT = new Key("dielectricConstant");
    //final Key EVALUATED_LOCATION_ENERGY = new Key("tetherEnergy");
    final Key TETHER_ENERGY = new Key("tetherEnergy");
    final Key CA_TETHER_ENERGY = new Key("caTetherEnergy");
    final Key CALPHA_HYDROGEN_BONDS = new Key("cAlphaHydrogenBonds");
    final Key CALPHA_HYDROGEN_BONDS_PLANE = new Key("cAlphaPlane");
    final Key HYDROGEN_BONDS_ANGLES = new Key("hydrogenBondsAngles");
    final Key HYDROGEN_BONDS_PLANE = new Key("hydrogenBondsPlane");
    final Key DISTANCE_CONSTRAINT_PCA = new Key("distanceConstraintPCA");
    final Key WARP_ENERGY = new Key("warpEnergy");
    final Key UN_WARP_ENERGY = new Key("unWarpEnergy");
    final Key LINEAR_RG = new Key("linearRG");
    final Key RG_OF_CHARGED_ATOMS = new Key("rgOfChargedAtoms");
    final Key RG = new Key("rg");
    final Key RG_LOGARITMIC = new Key("rgLogaritmicWeight");
    final Key RG_RATIO = new Key("rgRatioWeight");
    final Key RG_LINEAR_POLAR = new Key("rgPolarWeight");
    final Key RG_LINEAR_NON_POLAR = new Key("rgNonPolarWeight");
    final Key RG_LINEAR_BACKBONE = new Key("rgBackboneWeight");
    final Key SMOOTH_ROTAMER_LIBRARY_ENERGY = new Key("smoothRotamerLibrary");
    final Key RAMACHANDRAN_SIDECHAIN_ENERGY = new Key("ramachandranSidechain");
    final Key RAMACHANDRAN = new Key("ramachandran");
    final Key COOPERATIVE_RAMACHANDRAN_ENERGY = new Key("cooperativeRamachandranSidechain");
    final Key COOPERATIVE_RAMACHANDRAN_FILENAME = new Key("cooperativeRamachandranSidechainFile");
    final Key COOPERATIVE_Z_RAMACHANDRAN_ENERGY = new Key("cooperativeZRamachandranSidechain");
    final Key COOPERATIVE_Z_RAMACHANDRAN_FILENAME = new Key("cooperativeZRamachandranSidechainFile");
    final Key COOPERATIVE_Zstd_RAMACHANDRAN_ENERGY = new Key("cooperativeZstdRamachandranSidechain");

    final Key COMPOSITE_PROPENSITY_ENERGY = new Key("compositePropensity");
    final Key COOPERATIVE_PROPENSITY_ENERGY = new Key("cooperativePropensity");
    final Key COOPERATIVE_PROPENSITY_FILENAME = new Key("cooperativePropensityFile");
    final Key COOPERATIVE_Z_PROPENSITY_ENERGY = new Key("cooperativeZPropensity");
    final Key COOPERATIVE_Z_PROPENSITY_FILENAME = new Key("cooperativeZPropensityFile");
    final Key COOPERATIVE_Zstd_PROPENSITY_ENERGY = new Key("cooperativeZstdPropensity");
    final Key TEMPLATE_ENERGY = new Key("templateEnergy");
    final Key ATOMIC_PAIRWISE_PMF_SUMMA_ENERGY = new Key("atomicPairwisePMFSumma");
    final Key COOPERATIVE_ATOMIC_PAIRWISE_PMF_SUMMA_ENERGY = new Key("cooperativeAtomicPairwisePMFSumma");
    final Key COOPERATIVE_ATOMIC_PAIRWISE_PMF_SUMMA_FILENAME = new Key("cooperativeAtomicPairwisePMFSummaFile");
    final Key COOPERATIVE_PERATOM_SUMMA_FILENAME = new Key("cooperativePerAtomSummaFile");
    final Key COOPERATIVE_PERATOM_SUMMA_ENERGY = new Key("cooperativePerAtomSumma");
    final Key COOPERATIVE_Z_SUMMA_FILENAME = new Key("cooperativeZSummaFile");
    final Key COOPERATIVE_Z_SUMMA_ENERGY_NON_POLAR = new Key("cooperativeZSummaNonPolar");
    final Key COOPERATIVE_Z_SUMMA_ENERGY_POLAR = new Key("cooperativeZSummaPolar");
    final Key COOPERATIVE_Z_SUMMA_ENERGY_BACKBONE_OO_NN = new Key("cooperativeZSummaBackboneNNOO");
    final Key COOPERATIVE_Z_SUMMA_ENERGY_BACKBONE = new Key("cooperativeZSummaBackbone");
    final Key COOPERATIVE_Z_SUMMA_ENERGY_NEUTRAL = new Key("cooperativeZSummaNeutral");
    final Key COOPERATIVE_Z_SUMMA_ENERGY = new Key("cooperativeZSumma");
    final Key COOPERATIVE_Zstd_SUMMA_ENERGY = new Key("cooperativeZstdSumma");
    final Key COOPERATIVE_Z_SOLVATE_FILENAME = new Key("cooperativeZSolvateFile");
    final Key COOPERATIVE_Z_SOLVATE_ENERGY = new Key("cooperativeZSolvate");


    //------------------------ inflate --------------------------
    final Key RMS_TARGET = new Key("RmsTarget");
    //------------------------ template based distance constraints --------------------------   
    final Key INTRA_SEGMENT_FACTOR = new Key("intraSegmentFactor");
    final Key INTRA_SEGMENT_TOLERANCE = new Key("intraSegmentTolerance");
    final Key INTER_SEGMENT_FACTOR = new Key("interSegmentFactor");
    final Key INTER_SEGMENT_TOLERANCE = new Key("interSegmentTolerance");
    final Key SATURATION = new Key("saturation");
    final Key UNSATISFIED_CUTTOF = new Key("unsatisfiedCutoff");
    final Key UP_TO_CUTOFF = new Key("upToCutoff");
    final Key DISTANCE_CONSTRAINTS_MASK = new Key("constraint");

    //--------------------------------- Superimpose -------------------------------------
    final Key SUPERIMPOSE = new Key("superimpose");
    final Key REFERENCE = new Key("reference");
    final Key RMS_MODE = new Key("mode");
    final Key ALL_CA = new Key("allCa");


    //--------------------------------- minimization loop -------------------------------------
    final Key MINIMIZATION_LOOP = new Key("minimizationLoop");
    final Key ITERATIONS_CA = new Key("nIterationsCA");
    final Key ITERATIONS_BACKBONE = new Key("nIterationsBackbone");
    final Key ITERATIONS_ALLATOM = new Key("nIterationsAllAtoms");
    //--------------------------------- MCM -------------------------------------
    final Key MC_MINIMIZATION = new Key("MCM");
    final Key INITIAL_TEMPERATURE = new Key("initialTemperature");
    final Key FINAL_TEMPERATURE = new Key("finalTemperature");
    final Key MCM_PERTURBATION = new Key("mcmPerturbation");

    //--------------------------------- Homology Modeling -------------------------------------
    final Key TARGET_NAME = new Key("targetName");
    final Key TEMPLATE_NAME = new Key("templateName");
    final Key TARGET_FILE_PATH = new Key("targetFilePath");
    final Key ALINMENT_FILE_PATH = new Key("alinmentFilePath");
    final Key TEMPLATE_FILE_PATH = new Key("templateFilePath");
    final Key TEMPLATE_STRUCTURE = new Key("templateStructure");
    final Key TEMPLATE_DSSP = new Key("templateDssp");
    final Key TEMPLATE_TARGET_ALIGNMENT = new Key("templateTargetAlignment");
    final Key OUTPUT_FILE_PATH = new Key("outputFilePath");
    final Key OUTPUT_FILE_NAME = new Key("outputFileName");
    final Key SS_NAME = new Key("ssName");
    final Key LOOSEN_EDGE_LENGTH = new Key("loosenEdgeLength");
    final Key NON_FROZEN_BOND_DEPTH = new Key("nonFrozenBondDepth");
    final Key NON_FROZEN_RADIUS = new Key("nonFrozenRadius");
    final Key NUMBER_OF_MODELS = new Key("numberOfModels");
    final Key TARGET_SEQUENCE = new Key("targetSequence");
    final Key MAX_RUN_TIME = new Key("maxRunTime");
    final Key METHOD = new Key("method");
    final Key CASP_GROUP = new Key("caspGroup");
    final Key ADD_LOOPS = new Key("addLoops");
    final Key HHR_FILE_PATH = new Key("HHRFilePath");
    final Key RSR_FILE_PATH = new Key("RSRFile");


    //--------------------------------- Beautify --------------------------------
    final Key SHOTGUN_MODEL = new Key("shotgunModel");
    final Key CA_MODEL = new Key("caModel");
    final Key NUMBER_OF_CA_ITERATIONS = new Key("numberOfCaIterations");
    final Key WARP_THRESHOLD = new Key("warpThreshold");
    final Key WARP_STEP_SIZE = new Key("warpStepSize");
    final Key FREE_FINAL_MINIMIZATION = new Key("freeFinalMinimization");
    final Key NUMBER_OF_RUNS = new Key("numberOfRuns");
    final Key CA_CLASH_DISTANCE = new Key("caClashDistance");
    final Key CA_SHORT_DISTANCE = new Key("caShortDistance");
    final Key CA_LONG_DISTANCE = new Key("caLongDistance");
    final Key BEAUTIFY_PROBLEMATIC_RANGE = new Key("problematicRange");
    //--------------------------------- Refinement  --------------------------------
    final Key MODEL = new Key("model");
    final Key MODEL_DSSP = new Key("modelDssp");

    //--------------------------------- analysis ---------------------------------
    final Key DICTIONARY_KEY = new Key("COMMENT");
    final Key MESHILOG_KEY = new Key("MESHILOG");
    final Key KEY_KEY = new Key("T1");
    final Key VALUE_KEY = new Key("V1");

    //--------------------------------- Sequence -------------------------------------
    final Key AA_SEQUENCE = new Key("aa sequence");
    final Key SS_SEQUENCE = new Key("ss sequence");
    final Key ACCESIBILITY_SEQUENCE = new Key("accesibility sequence");
    //
    //--------------------------------- Helium Cluster Minimization  ------------------------
    final Key N_ATOMS = new Key("nAtoms");
    //--------------------------------- Distance Matrix -----------------------------
    final Key R_MAX = new Key("rMax");
    final Key BUFFER_SIZE = new Key("bufferSize");
    final Key GRID_EDGE = new Key("gridEdge");
    // -----------------------------  meshi.surface ------------------------------------
    final Key PDB_FILE = new Key("pdbFile");
    final Key KOEHL_FILE = new Key("koehlFile");
    //--------------------------------- Misc -------------------------------------
    final Key ON = new Key("isOn");
    final Key OFF = new Key("off");
    final Key END = new Key("end");
    final Key WEIGHT = new Key("weight");
    final Key INPUT_FILE = new Key("inputFile");
    final Key CUTOFF = new Key("cutoff");
    final Key NONE = new Key("none");
    final Key USE_FAST_ARCCOS = new Key("useFastArcCos");
    final Key DRESSER_FRAGMENTS = new Key("dresserFragments");
    final Key ROTAMER_LIBRARY = new Key("rotamerLibrary");
    final Key FIX_N_TERMINAL = new Key("fixNterminal");
    final Key FIX_C_TERMINAL = new Key("fixCterminal");
    final Key SEED = new Key("seed");
    final Key CORPUS_FILE_NAME = new Key("corpusFileName");
    final Key SCORE_WEIGHT = new Key("scoreWeight");
    final Key SELECTION_SCORE_WEIGHT = new Key("selectionScoreWeight");

    //------------------------- SymmetryComplex -----------------------
    final Key SYMMETRY_ENERGY = new Key("symmetryEnergy");
    final Key CYLINDER_ENERGY = new Key("cylinderEnergy");
    final Key EDM_ENERGY = new Key("EDMEnergy");
    final Key EDM_ENERGY_FILE_NAME = new Key("EDMEnergyFileName");
    final Key TOPOLOGY_MAP = new Key("topologyMap");
    final Key NUMBER_OF_CHAINS = new Key("numberOfChains");
    final Key NUMBER = new Key("number");
    final Key LOOP1 = new Key("loop1");
    final Key LOOP2 = new Key("loop2");
    final Key ANGLE_X = new Key("angleX");
    final Key ANGLE_Z = new Key("angleZ");
    final Key WIDTH_OF_HAIRPIN = new Key("theMostApplicableWidthOfHairpin");
    final Key MIN_WIDTH_OF_HAIRPIN = new Key("minWidthOfHairpin");
    final Key MAX_WIDTH_OF_HAIRPIN = new Key("maxWidthOfHairpin");
    final Key STEPS = new Key("accordion-like_pattern");
    final Key CONSTRICT = new Key("toPickTheMostNarrowLoopFromSuitable");
    final Key STRICT_CLASHES = new Key("strictClashes");
    final Key CHECK_INTERLOOP_DISTANCE = new Key("checkDistanceBetweenLoops");
}
