/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.InfoType;

/**

 */

public class
        SolvationInfoElement extends EnergyInfoElement {
    public final EnergyInfoElement scPolarInfo, bbPolarInfo, scCarbonInfo, bbCarbonInfo, stdInfo, entropyInfo, hbInfo, buriedHbinfo;
    //For statistics
    public final EnergyInfoElement bbPolarN, bbCarbonN, scPolarN, scCarbonN;

    public SolvationInfoElement(InfoType type, String comment, double weightSCPolarSolvate,
                                double weightSCCarbonSolvate,
                                double weightBBPolarSolvate,
                                double weightBBCarbonSolvate,
                                double weightStdSolvate) {
        super(type, comment, 1);
        getChildren().add(scPolarInfo  = new EnergyInfoElement(InfoType.SOLVATION_SC_POLAR, "Solvation term for sidechain polar getAtoms", weightSCPolarSolvate));
        getChildren().add(scCarbonInfo = new EnergyInfoElement(InfoType.SOLVATION_SC_CARBON, "Solvation term for sidechain carbon getAtoms", weightSCCarbonSolvate));
        getChildren().add(bbPolarInfo  = new EnergyInfoElement(InfoType.SOLVATION_BB_POLAR, "Solvation term for backbone  polar getAtoms", weightBBPolarSolvate));
        getChildren().add(bbCarbonInfo = new EnergyInfoElement(InfoType.SOLVATION_BB_CARBON, "Solvation term for backbone  carbon getAtoms in", weightBBCarbonSolvate));
        getChildren().add(hbInfo       = new EnergyInfoElement(InfoType.SOLVATION_HB, "Hydrogen bonds regardless of environment", 0.0001));
        getChildren().add(buriedHbinfo = new EnergyInfoElement(InfoType.SOLVATION_BURIED_HB, "Buried hydrogen bonds ", 0.0001));
        getChildren().add(stdInfo      = new EnergyInfoElement(InfoType.SOLVATION_STD, "Distribution of solvation energy over the protein", weightStdSolvate));
        getChildren().add(entropyInfo  = new EnergyInfoElement(InfoType.SOLVATION_ENTROPY, "Distribution of solvation energy over the protein", weightStdSolvate));
    //For statistics
        getChildren().add(bbPolarN     = new EnergyInfoElement(InfoType.BB_POLAR_N, "bbPolarN", 0.0));
        getChildren().add(bbCarbonN    = new EnergyInfoElement(InfoType.BB_CARBON_N, "bbPolarN", 0.0));
        getChildren().add(scPolarN     = new EnergyInfoElement(InfoType.SC_POLAR_N, "bbPolarN", 0.0));
        getChildren().add(scCarbonN    = new EnergyInfoElement(InfoType.SC_CARBON_N, "bbPolarN", 0.0));
    }
}
