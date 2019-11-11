/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.string.StringList;
import meshi.util.file.MeshiLineReader;
import meshi.parameters.MeshiPotential;
import meshi.parameters.ResidueType;
import meshi.parameters.AtomType;

import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 09/06/2009
 * Time: 10:08:45
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZPropensityParameters implements KeyWords, MeshiPotential {
    protected final double[] mean = new double[ResidueType.values().length];
    protected final double[] std = new double[ResidueType.values().length];

    public CooperativeZPropensityParameters(CommandList commands) {
        String parametersFileName = commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/meshiPotential/" + commands.firstWord(COOPERATIVE_Z_PROPENSITY_FILENAME).secondWord();
        Utils.println("CooperativeZPropensityParametersFile: " + parametersFileName);
        MeshiLineReader mlrParameters = new MeshiLineReader(parametersFileName);

        StringList lines = new StringList(mlrParameters);
        ResidueType type;
        for (String str : lines) {
            StringTokenizer line = new StringTokenizer(str);
            //line.nextToken(); // ignore the first word in the line.
            String name = line.nextToken();
            type = ResidueType.type(name);
            //minAvgSigma[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
            mean[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
            std[type.ordinal()] = (new Double(line.nextToken())).doubleValue();
        }
    }

}
