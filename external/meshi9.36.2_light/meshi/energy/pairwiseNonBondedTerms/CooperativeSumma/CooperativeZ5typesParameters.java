/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.string.StringList;
import meshi.util.file.MeshiLineReader;
import meshi.parameters.MeshiPotential;
import meshi.parameters.AtomType;

import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 22/12/2009
 * Time: 13:17:51
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZ5typesParameters implements KeyWords, MeshiPotential {
    final static int NumOfClusters = 3;
    protected final double[][] mean = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] std = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] proportion = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] meanPolar = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] stdPolar = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] proportionPolar = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] meanNonPolar = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] stdNonPolar = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] proportionNonPolar = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] meanNeutral = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] stdNeutral = new double[NumOfClusters][AtomType.values().length];
    protected final double[][] proportionNeutral = new double[NumOfClusters][AtomType.values().length];
    public final double[][] meanNonPolarAndNeutral = new double[NumOfClusters][AtomType.values().length];
    public final double[][] meanPOLARS_BB = new double[NumOfClusters][AtomType.values().length];
    public final double[][] meanPOLARSnoH_NNorOO = new double[NumOfClusters][AtomType.values().length];
    public final double[][] stdNonPolarAndNeutral = new double[NumOfClusters][AtomType.values().length];
    public final double[][] stdPOLARS_BB = new double[NumOfClusters][AtomType.values().length];
    public final double[][] stdPOLARSnoH_NNorOO = new double[NumOfClusters][AtomType.values().length];
    public final double[][] proportionNonPolarAndNeutral = new double[NumOfClusters][AtomType.values().length];
    public final double[][] proportionPOLARS_BB = new double[NumOfClusters][AtomType.values().length];
    public final double[][] proportionPOLARSnoH_NNorOO = new double[NumOfClusters][AtomType.values().length];


    public CooperativeZ5typesParameters(CommandList commands) {
        String parametersFileName = commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/meshiPotential/" + commands.firstWord(COOPERATIVE_Z_SUMMA_FILENAME).secondWord();

        Utils.println("CooperativeZSummaParametersFile: " + parametersFileName);
        MeshiLineReader mlrParameters = new MeshiLineReader(parametersFileName);

        StringList lines = new StringList(mlrParameters);
        AtomType type;
        for (String str : lines) {
            StringTokenizer line = new StringTokenizer(str);
            //line.nextToken(); // ignore the first word in the line.
            String name = line.nextToken();
            type = AtomType.type(name);
            for (int i = 0; i < NumOfClusters; i++) {
                mean[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                std[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                proportion[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
//          if (std[i][type.ordinal()] != 0)
                //          proportion[i][type.ordinal()] = proportion[i][type.ordinal()]/std[i][type.ordinal()];
            }
            for (int i = 0; i < NumOfClusters; i++) {
                meanPolar[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                stdPolar[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                proportionPolar[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                //        if (stdPolar[i][type.ordinal()] != 0)
                //        proportionPolar[i][type.ordinal()] = proportionPolar[i][type.ordinal()]/stdPolar[i][type.ordinal()];
            }
            for (int i = 0; i < NumOfClusters; i++) {
                meanNonPolar[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                stdNonPolar[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                proportionNonPolar[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                //    if (stdNonPolar[i][type.ordinal()] != 0)
                //     proportionNonPolar[i][type.ordinal()] = proportionNonPolar[i][type.ordinal()]/stdNonPolar[i][type.ordinal()];
            }
            for (int i = 0; i < NumOfClusters; i++) {
                meanNeutral[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                stdNeutral[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                proportionNeutral[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
//            if (stdNeutral[i][type.ordinal()] != 0)
                //             proportionNeutral[i][type.ordinal()] = proportionNeutral[i][type.ordinal()]/stdNeutral[i][type.ordinal()];
            }
            for (int i = 0; i < NumOfClusters; i++) {
                meanNonPolarAndNeutral[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                stdNonPolarAndNeutral[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                proportionNonPolarAndNeutral[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
//            if (stdNonPolarAndNeutral[i][type.ordinal()] != 0)
                //             proportionNonPolarAndNeutral[i][type.ordinal()] = proportionNonPolarAndNeutral[i][type.ordinal()]/stdNonPolarAndNeutral[i][type.ordinal()];
            }
            for (int i = 0; i < NumOfClusters; i++) {
                meanPOLARSnoH_NNorOO[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                stdPOLARSnoH_NNorOO[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                proportionPOLARSnoH_NNorOO[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
//            if (stdPOLARSnoH_NNorOO[i][type.ordinal()] != 0)
                //             proportionPOLARSnoH_NNorOO[i][type.ordinal()] = proportionPOLARSnoH_NNorOO[i][type.ordinal()]/stdPOLARSnoH_NNorOO[i][type.ordinal()];
            }
            for (int i = 0; i < NumOfClusters; i++) {
                meanPOLARS_BB[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                stdPOLARS_BB[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
                proportionPOLARS_BB[i][type.ordinal()] = (new Double(line.nextToken())).doubleValue();
//            if (stdPOLARS_BB[i][type.ordinal()] != 0)
                //             proportionPOLARS_BB[i][type.ordinal()] = proportionPOLARS_BB[i][type.ordinal()]/stdPOLARS_BB[i][type.ordinal()];
            }
        }
    }

}

