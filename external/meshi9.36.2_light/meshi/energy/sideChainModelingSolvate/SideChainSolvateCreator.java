/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.sideChainModelingSolvate;

import meshi.energy.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.info.InfoType;
import meshi.util.string.*;

import java.util.*;

import meshi.util.*;

/**
 * This class is a truncated form of the corresponding class in the "meshi.energy.solvation" package. It was designed
 * solely for accelerating the application SCMOD (concurrent sidechain modeling), and should not be used for other
 * purposes.
 * <p/>
 * Please do not use it!
 */


public class SideChainSolvateCreator extends EnergyCreator implements KeyWords {
    private static ArrayList parametersList = null;
    // Default values for the HB energy term weight and the HB angular score thresholds.  
    private double simpleHBweight = 0.0;
    private double sigmoidBeginsWithH = 90.0;
    private double sigmoidEndsWithH = 100.0;
    private double sigmoidBeginsNoH = 80.0;
    private double sigmoidEndsNoH = 90.0;

    public SideChainSolvateCreator(double cooperativeSolvateWeight,
                                   double simpleHBweight,
                                   double sigmoidBeginsWithH,
                                   double sigmoidEndsWithH,
                                   double sigmoidBeginsNoH,
                                   double sigmoidEndsNoH) {
        this(cooperativeSolvateWeight, simpleHBweight);
        this.sigmoidBeginsWithH = sigmoidBeginsWithH;
        this.sigmoidEndsWithH = sigmoidEndsWithH;
        this.sigmoidBeginsNoH = sigmoidBeginsNoH;
        this.sigmoidEndsNoH = sigmoidEndsNoH;
    }

    public SideChainSolvateCreator(double cooperativeSolvateWeight, double simpleHBweight) {
        this(cooperativeSolvateWeight);
        this.simpleHBweight = simpleHBweight;
    }

    public SideChainSolvateCreator(double weight) {
        super(weight,InfoType.SOLVATION);
    }

    public SideChainSolvateCreator() {
        super(InfoType.SOLVATION);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if (parametersList == null) {
            // Appanding the path to the list of parameter filenames.
            String[] strlist = new String[SOLVATE_SC_PARAMETERS.length];
            String pathname = parametersDirectory(commands).concat("/");
            for (int cc = 0; cc < SOLVATE_SC_PARAMETERS.length; cc++)
                strlist[cc] = pathname.concat(SOLVATE_SC_PARAMETERS[cc]);
            parametersList = new SideChainSolvateParametersList(strlist);
        }
        EnergyInfoElement info = new EnergyInfoElement(InfoType.SOLVATION, " Solvation energy for side chain modeling", weight);
        SideChainSolvateEnergy solvateEnergy = new SideChainSolvateEnergy(protein.atoms(),
                distanceMatrix,
                (SideChainSolvateParametersList) parametersList,
                simpleHBweight,
                sigmoidBeginsWithH,
                sigmoidEndsWithH,
                sigmoidBeginsNoH,
                sigmoidEndsNoH,
                info);
        solvateEnergy.setComment("Solvate");
        return solvateEnergy;
    }

}
