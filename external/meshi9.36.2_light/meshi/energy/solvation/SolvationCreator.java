/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * The different between this creator and the "SolvateCreatorRegularHB" (see the latter main comment for a lot of details), is in the
 * way the hydrogen-bond distance sigmoid is setResidue. Here parameters are chosen so that the hydrogen-bond starts its drop to 0 at a
 * donor-acceptor distance of 3.1 Ang, and is 0 only for donor-acceptor distance of 4.0 Ang. This is good to increase the basin of
 * attraction of direct minimizations protocols.
 */

public class SolvationCreator extends EnergyCreator {

    // The different weights relevent to this term.
    protected double weightSCPolarSolvate = -9999.0;
    protected double weightSCCarbonSolvate = -9999.0;
    protected double weightBBPolarSolvate = -9999.0;
    protected double weightBBCarbonSolvation = -9999.0;
    protected double weightStdSolvation = -9999.0;
    private static SolvateParametersList parametersList = null;
    private String comment = null;
    private DistanceMatrix distanceMatrix = null;

    public SolvationCreator(double weightSCPolarSolvate, double weightSCCarbonSolvate,
                            double weightBBPolarSolvate, double weightBBCarbonSolvation, double weightHB) {
        this(weightSCPolarSolvate, weightSCCarbonSolvate, weightBBPolarSolvate, weightBBCarbonSolvation, weightHB, "Solvation");
    }

    public SolvationCreator(double weight, double weightSCPolarSolvate, double weightSCCarbonSolvate,
                            double weightBBPolarSolvate, double weightBBCarbonSolvation, double weightHB) {
        this(weightSCPolarSolvate, weightSCCarbonSolvate, weightBBPolarSolvate, weightBBCarbonSolvation, weightHB, "Solvation");
        this.weight = weight;
    }

    public SolvationCreator(double weightSCPolarSolvate, double weightSCCarbonSolvate,
                            double weightBBPolarSolvate, double weightBBCarbonSolvation, double weightStdSolvation, String comment) {
        super(InfoType.SOLVATION);
        this.weightSCPolarSolvate = weightSCPolarSolvate;
        this.weightSCCarbonSolvate = weightSCCarbonSolvate;
        this.weightBBPolarSolvate = weightBBPolarSolvate;
        this.weightBBCarbonSolvation = weightBBCarbonSolvation;
        this.weightStdSolvation = weightStdSolvation;
        this.comment = comment;
    }

    public SolvationCreator(double weight) {
        super(weight,InfoType.SOLVATION);
        this.weightSCPolarSolvate = weight;
        this.weightSCCarbonSolvate = weight;
        this.weightBBPolarSolvate = weight;
        this.weightBBCarbonSolvation = weight;
        this.weightStdSolvation = weight;
    }

    public SolvationCreator() {
        super(InfoType.SOLVATION);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if ((term != null) && (this.distanceMatrix == distanceMatrix)) //chen 3/3/12
                return term;
            this.distanceMatrix = distanceMatrix;
            if (weightSCPolarSolvate < 0) {
                weightSCPolarSolvate        = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_SC_POLAR.tag).thirdWordDouble();
                weightSCCarbonSolvate       = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_SC_CARBON.tag).thirdWordDouble();
                weightBBPolarSolvate        = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_BB_POLAR.tag).thirdWordDouble();
                weightBBCarbonSolvation     = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_BB_CARBON.tag).thirdWordDouble();
                weightStdSolvation          = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_STD.tag).thirdWordDouble();
                weightSCPolarSolvate       *= weight();
                weightSCCarbonSolvate      *= weight();
                weightBBPolarSolvate       *= weight();
                weightBBCarbonSolvation    *= weight();
                weightStdSolvation         *= weight();
            }
        if (parametersList == null) {
            // Appanding the path to the list of parameter filenames.
            String[] strlist = new String[SOLVATION_MINIMIZE_HB_PARAMETERS.length];
            String pathname = parametersDirectory(commands).concat("/");
            for (int cc = 0; cc < SOLVATION_MINIMIZE_HB_PARAMETERS.length; cc++)
                strlist[cc] = pathname.concat(SOLVATION_MINIMIZE_HB_PARAMETERS[cc]);
            parametersList = new SolvateParametersList(strlist);
        }



        return term = createTerm(protein.atoms(),
                distanceMatrix,
                parametersList);
    }

    protected AbstractEnergy createTerm(AtomList atomList, DistanceMatrix distanceMatrix, SolvateParametersList parametersList) {
        SolvationInfoElement info = new SolvationInfoElement(InfoType.SOLVATION, "Solvation energy",
                weightSCPolarSolvate,
                weightSCCarbonSolvate,
                weightBBPolarSolvate,
                weightBBCarbonSolvation,
                weightStdSolvation);
        return new SolvationEnergy(atomList,
                distanceMatrix,
                parametersList,
                info);
    }

}
