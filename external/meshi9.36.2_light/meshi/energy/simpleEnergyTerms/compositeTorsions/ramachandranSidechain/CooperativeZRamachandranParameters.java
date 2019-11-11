/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.string.StringList;
import meshi.util.file.MeshiLineReader;
import meshi.parameters.MeshiPotential;
import meshi.parameters.ResidueType;

import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 09/06/2009
 * Time: 13:20:56
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZRamachandranParameters implements KeyWords, MeshiPotential {
    protected final double[] mean = new double[ResidueType.values().length];
    protected final double[] std = new double[ResidueType.values().length];


    public CooperativeZRamachandranParameters(CommandList commands) {
        String parametersFileName = commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/meshiPotential/" + commands.firstWord(COOPERATIVE_Z_RAMACHANDRAN_FILENAME).secondWord();
        Utils.println("CooperativeZRamachandranParametersFile: " + parametersFileName);
        MeshiLineReader mlrParameters = new MeshiLineReader(parametersFileName);

        StringList lines = new StringList(mlrParameters);
        ResidueType type;
        for (String str : lines) {
            StringTokenizer line = new StringTokenizer(str);
            String name = line.nextToken();
            type = ResidueType.type(name);
            mean[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
            std[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
        }
    }
}

