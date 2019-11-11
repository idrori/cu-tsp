/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.parameters.ResidueType;
import meshi.parameters.MeshiPotential;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.string.StringList;
import meshi.util.file.MeshiLineReader;

import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 09/02/2009
 * Time: 10:33:34
 * To change this template use File | Settings | File Templates.
 */
public class CooperativePropensityParameters implements KeyWords, MeshiPotential {
    protected final double[] minAvgSigma = new double[ResidueType.values().length];

    public CooperativePropensityParameters(CommandList commands) {
        String parametersFileName = commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/meshiPotential/" + commands.firstWord(COOPERATIVE_PROPENSITY_FILENAME).secondWord();
        System.out.println("CooperativePropensityParametersFile: " + parametersFileName);
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
