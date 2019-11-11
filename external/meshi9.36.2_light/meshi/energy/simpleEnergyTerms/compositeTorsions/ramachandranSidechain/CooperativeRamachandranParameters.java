/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.simpleEnergyTerms.compositeTorsions.*;

import meshi.parameters.*;
import meshi.util.string.*;
import meshi.util.*;
import meshi.util.file.*;

import java.util.*;

public class CooperativeRamachandranParameters implements KeyWords, MeshiPotential {
    protected final double[] minAvgSigma = new double[ResidueType.values().length];

    public CooperativeRamachandranParameters(CommandList commands) {
        //String          parametersFileName = commands.firstWord( PARAMETERS_DIRECTORY ).secondWord() + "/" + COOPERATIVE_RAMACHANDRAN_PARAMETERS;
        String parametersFileName = commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/meshiPotential/" + commands.firstWord(COOPERATIVE_RAMACHANDRAN_FILENAME).secondWord();
        System.out.println("CooperativeRamachandranParametersFile: " + parametersFileName);
        MeshiLineReader mlrParameters = new MeshiLineReader(parametersFileName);
        StringList lines = new StringList(mlrParameters);
        ResidueType type;
        for (String str : lines) {
            StringTokenizer line = new StringTokenizer(str);
            line.nextToken(); // ignore the first word in the line.
            String name = line.nextToken();
            type = ResidueType.type(name);
            minAvgSigma[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
        }
    }
}