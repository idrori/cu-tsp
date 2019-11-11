/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.parameters;

public interface MeshiPotential {
    public final String HYDROGEN_BONDS_PAIRS_PARAMETERS_SURFACE = "meshiPotential/HydrogenBondsParameters.meshi.surface.dat";
    public final String HYDROGEN_BONDS_PAIRS_BETA_PARAMETERS = "meshiPotential/HydrogenBondsBetaParameters.dat";
    public final String HYDROGEN_BONDS_PAIRS_HELIX_PARAMETERS = "meshiPotential/HydrogenBondsHelixParameters.dat";
    public final String HYDROGEN_BONDS_PAIRS_BETA_BACKGROUND_PARAMETERS = "meshiPotential/HydrogenBondsBetaBackgroundParameters.dat";
    public final String HYDROGEN_BONDS_PAIRS_HELIX_BACKGROUND_PARAMETERS = "meshiPotential/HydrogenBondsHelixBackgroundParameters.dat";
    public final String BOND_PARAMETERS = "meshiPotential/bondEnergyParameters.dat";
    public final String ANGLE_PARAMETERS = "meshiPotential/angleEnergyParameters.dat";
    public final String PLANE_PARAMETERS = "meshiPotential/planeEnergyParameters.dat";
    public final String OUT_OF_PLANE_PARAMETERS = "meshiPotential/outOfPlaneEnergyParameters.dat";
    public final String LENNARD_JONES_PARAMETERS = "meshiPotential/LennardJonesParameters.dat";
    public final String LENNARD_JONES_PARAMETERS_CA = "meshiPotential/LennardJonesParametersCa.dat";
    public final String LENNARD_JONES_PARAMETERS_BACKBONE = "meshiPotential/LennardJonesParametersBackbone.dat";
    public final String CONTACTS_PARAMETERS = "meshiPotential/contactsParameters.dat";
    public final String LJ_ENVIRONMENT_PARAMETERS = "meshiPotential/LJEnvironment.dat";
    public final String LJ_ENVIRONMENT_PARAMETERS_CA = "meshiPotential/LJEnvironmentCa.dat";
    public final String LJ_ENVIRONMENT_PARAMETERS_BACKBONE = "meshiPotential/LJEnvironmentBackbone.dat";
    public final String CONTACTS_ENVIRONMENT_PARAMETERS = "meshiPotential/contactsEnvironment.dat";
    public final String EXCLUDED_VOL_PARAMETERS = "meshiPotential/ExcludedVolumeParameters.dat";

    /* parameters file name for the atomicPairwisePMFSumma energy term */
    public final String ATOMIC_PAIRWISE_PMF_SUMMA_PARAMETERS = "meshiPotential/PairwiseAtomicPMF_summa_5.5.dat";
    public final String COOPERATIVE_ATOMIC_PAIRWISE_PMF_SUMMA_PARAMETERS = "meshiPotential/PairwiseAtomicPMF_summa_5.5.cooperative.dat";
    // public final String COOPERATIVE_ATOMIC_PAIRWISE_PMF_SUMMA_PARAMETERS = "meshiPotential/cooperative_summa_1_without0.dat";
    /* parameters file name for the compositeTorsionsEnergy term */
    public final String COMPOSITE_TORSIONS_PARAMETERS = "meshiPotential/CompositeTorsionsParameters.dat";
    public final String COMPOSITE_PROPENSITY_PARAMETERS = "meshiPotential/CompositePropensity2D.dat";
    public final String COMPOSITE_PROPENSITY_2D_PARAMETERS = "meshiPotential/CompositePropensity2D_withPP.dat";

    public final String ELECTROSTATICS_PARAMETERS = "meshiPotential/ElectrostaticParameters.dat";
    public final String[] TWO_TORSIONS_PARAMETERS = {
            "meshiPotential/parametersPhiChiCoil.dat",
            "meshiPotential/parametersPhiChiHelix.dat",
            "meshiPotential/parametersPhiChiStrand.dat",

            "meshiPotential/parametersPsiChiCoil.dat",
            "meshiPotential/parametersPsiChiHelix.dat",
            "meshiPotential/parametersPsiChiStrand.dat",


            "meshiPotential/parametersCHI1CHI2.dat",
            "meshiPotential/parametersCHI2CHI3.dat",
            "meshiPotential/parametersCHI3CHI4.dat",

            "meshiPotential/parametersPhiPsiCoil.dat",
            "meshiPotential/parametersPhiPsiHelix.dat",
            "meshiPotential/parametersPhiPsiStrand.dat"};

    public final String COMPOSITE_PROPENSITY_2D_WITH_PP_PARAMETERS = "meshiPotential/CompositePropensity2D_withPP.dat";
    //public final String COOPERATIVE_PROPENSITY_PARAMETERS = "meshiPotential/cooperativePropensity_without0123.dat";
    public final String COOPERATIVE_PROPENSITY_PARAMETERS = "meshiPotential/cooperativePropensity.dat";
    //public final String COOPERATIVE_PROPENSITY_PARAMETERS = "meshiPotential/cooperative_propen_1_without0.dat";
    public final String COOPERATIVE_RAMACHANDRAN_PARAMETERS = "meshiPotential/CooperativeRamachandran.dat";
    //public final String COOPERATIVE_RAMACHANDRAN_PARAMETERS = "meshiPotential/cooperative_ramach_1_without0.dat";
    public final String[] FLAT_RAMACH_PARAMETERS = {
            "meshiPotential/flatRamachHelix.dat",
            "meshiPotential/flatRamachStrand.dat"};
    public final String[] FLAT_RAMACH_PARAMETERS_COIL = {
            "meshiPotential/flatRamachCoil.dat"};

    public final String[] PROPENSITY_TORSION_PARAMETERS = {"meshiPotential/parametersPropensity.dat"};

    public final String[] PROPENSITY_ANGLE_PARAMETERS = {"meshiPotential/parametersPropensityAngle.dat"};

    public final String[] ONE_ANGLE_PARAMETERS = {"meshiPotential/parameters3CAAll.dat",
            "meshiPotential/parameters3CACoil.dat",
            "meshiPotential/parameters3CAHelix.dat",
            "meshiPotential/parameters3CAStrand.dat"};

    public final String[] TWO_ANGLES_PARAMETERS = {"meshiPotential/parameters3CACONSECAll.dat",
            "meshiPotential/parameters3CACONSECCoil.dat",
            "meshiPotential/parameters3CACONSECHelix.dat",
            "meshiPotential/parameters3CACONSECStrand.dat"};

    public final String ALPHA_ANGLE_PARAMETERS = "meshiPotential/alphaRangeParameters.dat";

    public final String ALPHA_TORSION_PARAMETERS = "meshiPotential/alphaTorsionRangeParameters.dat";

    public final String HELIX = "HELIX"; //Effects both TwoTiorsions and HydrogenBondsPairs
    public final String SHEET = "SHEET"; //Effects both TwoTiorsions and HydrogenBondsPairs
    public final String COIL = "COIL";  //Effects both TwoTiorsions and HydrogenBondsPairs
    public final String HELIX_OR_COIL = "HELIX_OR_COIL"; // Effects only HydrogenBondsPairs
    public final String SHEET_OR_COIL = "SHEET_OR_COIL"; // Effects only HydrogenBondsPairs

    public final String ACCESSIBLE = "ACCESSIBLE";
    public final String BURIED = "BURIED";


    public final String[] SOLVATE_PARAMETERS = {"meshiPotential/SolvateCend.dat",
            "meshiPotential/SolvateCp1.dat",
            "meshiPotential/SolvateCp2.dat",
            "meshiPotential/SolvateCvalAtp1.dat",
            "meshiPotential/SolvateCvalAtp2.dat",
            "meshiPotential/SolvateRegularHBend.dat",
            "meshiPotential/SolvateRegularHBp1.dat",
            "meshiPotential/SolvateRegularHBp2.dat",
            "meshiPotential/SolvateRegularHBvalAtp1.dat",
            "meshiPotential/SolvateRegularHBvalAtp2.dat",
            "meshiPotential/SolvateMESHI2Tsai.dat",
            "meshiPotential/SolvatePolarSideChainSplines.dat",
            "meshiPotential/SolvateCarbonSideChainSplines.dat",
            "meshiPotential/SolvatePolarBackboneSplines.dat"};

    public final String[] SOLVATE_LONG_HB_PARAMETERS = {"meshiPotential/SolvateCend.dat",
            "meshiPotential/SolvateCp1.dat",
            "meshiPotential/SolvateCp2.dat",
            "meshiPotential/SolvateCvalAtp1.dat",
            "meshiPotential/SolvateCvalAtp2.dat",
            "meshiPotential/SolvateLongHBend.dat",
            "meshiPotential/SolvateLongHBp1.dat",
            "meshiPotential/SolvateLongHBp2.dat",
            "meshiPotential/SolvateLongHBvalAtp1.dat",
            "meshiPotential/SolvateLongHBvalAtp2.dat",
            "meshiPotential/SolvateMESHI2Tsai.dat",
            "meshiPotential/SolvatePolarSideChainSplines.dat",
            "meshiPotential/SolvateCarbonSideChainSplines.dat",
            "meshiPotential/SolvatePolarBackboneSplines.dat"};

    public final String[] SOLVATE_MINIMIZE_HB_PARAMETERS = {"meshiPotential/SolvateCend.dat",
            "meshiPotential/SolvateCp1.dat",
            "meshiPotential/SolvateCp2.dat",
            "meshiPotential/SolvateCvalAtp1.dat",
            "meshiPotential/SolvateCvalAtp2.dat",
            "meshiPotential/SolvateForMinimizeHBend.dat",
            "meshiPotential/SolvateForMinimizeHBp1.dat",
            "meshiPotential/SolvateForMinimizeHBp2.dat",
            "meshiPotential/SolvateForMinimizeHBvalAtp1.dat",
            "meshiPotential/SolvateForMinimizeHBvalAtp2.dat",
            "meshiPotential/SolvateMESHI2Tsai.dat",
            "meshiPotential/SolvatePolarSideChainSplines.dat",
            "meshiPotential/SolvateCarbonSideChainSplines.dat",
            "meshiPotential/SolvatePolarBackboneSplines.dat"};

    public final String[] SOLVATION_MINIMIZE_HB_PARAMETERS = {"meshiPotential/SolvateCend.dat",
            "meshiPotential/SolvateCp1.dat",
            "meshiPotential/SolvateCp2.dat",
            "meshiPotential/SolvateCvalAtp1.dat",
            "meshiPotential/SolvateCvalAtp2.dat",
            "meshiPotential/SolvateForMinimizeHBend.dat",
            "meshiPotential/SolvateForMinimizeHBp1.dat",
            "meshiPotential/SolvateForMinimizeHBp2.dat",
            "meshiPotential/SolvateForMinimizeHBvalAtp1.dat",
            "meshiPotential/SolvateForMinimizeHBvalAtp2.dat",
            "meshiPotential/SolvateMESHI2Tsai.dat",
            "meshiPotential/SolvationSplineParameters.dat"};

    public final String[] SOLVATE_NOHB_PARAMETERS = {"meshiPotential/SolvateMESHI2Tsai.dat",
            "meshiPotential/SolvateNoHydrogenBondParameters.dat",
            "meshiPotential/SolvateExtResCend.dat",
            "meshiPotential/SolvateExtResCp1.dat",
            "meshiPotential/SolvateExtResCp2.dat",
            "meshiPotential/SolvateExtResCvalAtp1.dat",
            "meshiPotential/SolvateExtResCvalAtp2.dat"};
    public final String[] SOLVATE_SC_PARAMETERS = {"meshiPotential/sc/SolvateMESHI2Tsai.dat",
            "meshiPotential/sc/SolvateExtResParameters.dat",
            "meshiPotential/sc/SolvateExtResCend.dat",
            "meshiPotential/sc/SolvateExtResCp1.dat",
            "meshiPotential/sc/SolvateExtResCp2.dat",
            "meshiPotential/sc/SolvateExtResCvalAtp1.dat",
            "meshiPotential/sc/SolvateExtResCvalAtp2.dat",
            "meshiPotential/sc/SolvateExtResHBend.dat",
            "meshiPotential/sc/SolvateExtResHBp1.dat",
            "meshiPotential/sc/SolvateExtResHBp2.dat",
            "meshiPotential/sc/SolvateExtResHBvalAtp1.dat",
            "meshiPotential/sc/SolvateExtResHBvalAtp2.dat"};

}
