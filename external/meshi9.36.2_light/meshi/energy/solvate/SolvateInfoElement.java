/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

/**

 */

public class SolvateInfoElement extends EnergyInfoElement {
    public final MeshiInfo scPolarInfo, bbPolarInfo, scCarbonInfo, stdInfo;

    public SolvateInfoElement(InfoType type, String comment, double weightSCPolarSolvate,
                              double weightSCCarbonSolvate,
                              double weightBBPolarSolvate) {
        super(type, comment, 1);
        scPolarInfo  = new EnergyInfoElement(InfoType.SOLVATE_SC_POLAR, "Solvation term for sidechain polar atoms", weightSCPolarSolvate);
        scCarbonInfo = new EnergyInfoElement(InfoType.SOLVATE_SC_CARBON, "Solvation term for sidechain carbon atoms", weightSCCarbonSolvate);
        bbPolarInfo  = new EnergyInfoElement(InfoType.SOLVATE_BB_POLAR, "Solvation term for backbone  polar atoms", weightBBPolarSolvate);
        stdInfo      = new EnergyInfoElement(InfoType.SOLVATE_STD," Distribution of solvation energy over the protein");
        getChildren().add(scPolarInfo);
        getChildren().add(scCarbonInfo);
        getChildren().add(bbPolarInfo);
        getChildren().add(stdInfo);
    }

    MeshiInfo stdElememt() {
        return stdInfo;
    }
}
