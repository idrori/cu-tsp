/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.scoringFunctions;

import meshi.util.KeyWords;
import meshi.util.Utils;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 10/11/2010
 * Time: 23:18:00
 * To change this template use File | Settings | File Templates.
 */
public class CombinedEnergyScore implements Score, KeyWords {
    public final Comparator<NormalizedEnergyScore> weightComparator = new WeightComparator();
    private NormalizedEnergyScore[] scoreFunctions;
    private double interdecileRange = -9999;
//    private TotalEnergy energy;
    private double[] cumulativeSumWeights;
    private DoubleInfoElement lengthElement, changeElement;
    private String name;
    private String[] fields;

    public CombinedEnergyScore(String parametersFileName, String name) {
        this.name = name;
        File parametersFile = new File(parametersFileName);
        ConfigurationArray configurationArray= new ConfigurationArray(parametersFile);
        scoreFunctions = new NormalizedEnergyScore[configurationArray.size()];
        int iConfig = 0;
        fields = configurationArray.get(0).fields;
        for (Configuration configuration:configurationArray)     {
            scoreFunctions[iConfig] = new NormalizedEnergyScore(configuration);
            iConfig++;
        }
        cumulativeSumWeights = new double[configurationArray.size()];
        lengthElement = null;
        changeElement = null;
    }

    public String toString() {
        return name;
    }

    public void setLengthElement(int length){
        lengthElement = new DoubleInfoElement(InfoType.LENGTH,"The length of the protein",length);
    }
    public void setChangeElement(double change) {
        changeElement = new DoubleInfoElement(InfoType.CHANGE, "structural change (in RMS) due to refinement", change);
    }

    public boolean fieldRequired(String s) {
        for (String field :fields)
            if (field.equals(s)) return true;
        return false;
    }
    public MeshiInfo score(MeshiInfo energyInfo) {
        MeshiInfo out;

        for (int i = scoreFunctions.length - 1; i >= 0; i--) {// Makes it easier to compare to the MATLAB version
            NormalizedEnergyScore scoreFunction = scoreFunctions[i];
            if (scoreFunction == null)
                throw new RuntimeException("This is weird " + i);
            scoreFunction.calcScore(energyInfo.flatten());
        }
        if (scoreFunctions.length == 1) {
            interdecileRange = 0;
            out = new MeshiInfo(InfoType.SIMPLE_SCORE, scoreFunctions[0].getScore(), "A single predictor score.");
            return out;
        }

        Arrays.sort(scoreFunctions);
        double sumWeights = 0;
        for (int i = 0; i < cumulativeSumWeights.length; i++)
            cumulativeSumWeights[i] = sumWeights = sumWeights + scoreFunctions[i].getWeight();

        MeshiInfo[] scores = new MeshiInfo[3];
        out = new MeshiInfo(InfoType.WEIGHTED_MEDIAN_SCORE,
                getMedian(scoreFunctions, cumulativeSumWeights, sumWeights),
                toString()+"_"+InfoType.WEIGHTED_MEDIAN_SCORE.tag);
        int top10Index = Arrays.binarySearch(cumulativeSumWeights, 9 * sumWeights / 10);
        if (top10Index < 0) top10Index = -top10Index;
        int bottom10Index = Arrays.binarySearch(cumulativeSumWeights, sumWeights / 10);
        if (bottom10Index < 0) bottom10Index = -bottom10Index;
        double interdecile = scoreFunctions[top10Index].getScore()-scoreFunctions[bottom10Index].getScore();

        Arrays.sort(scoreFunctions,weightComparator);
        top10Index = 9*scoreFunctions.length/10;
        NormalizedEnergyScore[] topTenPredictors = new NormalizedEnergyScore[scoreFunctions.length - top10Index];
        double[] topTenSumWeigths = new double[scoreFunctions.length - top10Index];
        sumWeights = 0;
        for (int i = top10Index; i < scoreFunctions.length; i++) {
            topTenPredictors[i - top10Index] = scoreFunctions[i];
        }
        Arrays.sort(topTenPredictors);
        for (int i = 0; i < topTenSumWeigths.length; i++)
            topTenSumWeigths[i] = sumWeights = sumWeights + topTenPredictors[i].getWeight();

        scores[0] = new MeshiInfo(InfoType.TOP_TEN_AVERAGE_SCORE, getMean(topTenPredictors, sumWeights),
                                          "Average score of the top 10% predictors");
        scores[1] = new MeshiInfo(InfoType.TOP_TEN_MEDIAN_SCORE, getMedian(topTenPredictors, topTenSumWeigths, sumWeights),
                                          "Weighted median of the top 10% predictors' scores");
        scores[2] = new MeshiInfo(InfoType.INTERDECILE, interdecile,
                    toString()+"_"+InfoType.INTERDECILE.tag);
        out.getChildren().add(scores[2]);
        return out;
    }

    public static double getMean(NormalizedEnergyScore[] scoreFunctions, double sumWeights) {
        double sum = 0;
        for (NormalizedEnergyScore scoreFunction : scoreFunctions) {
            sum += scoreFunction.getScore()*scoreFunction.getWeight();
        }
        return sum/sumWeights;
    }
    public static double getMedian(NormalizedEnergyScore[] scoreFunctions, double[] cumulativeSumWeights,double sumWeights) {
        int midPoint = Arrays.binarySearch(cumulativeSumWeights, sumWeights / 2);
        if (midPoint < 0) //sumWeights/2 is not explicitly in the array and binarry search provide the "insertion point"
            midPoint = -midPoint;

        int topMedianIndex = midPoint; //It is the index above  sumWeights/2
        double topMedian = scoreFunctions[topMedianIndex].getScore();
        int bottomMedianIndex = midPoint - 1;
        double bottomWeight = cumulativeSumWeights[topMedianIndex] - sumWeights / 2;
        double topWeight;
        double bottomMedian;
        if (topMedianIndex == 0) {
            bottomMedian = 0;
            topWeight = cumulativeSumWeights[topMedianIndex];
        } else {
            bottomMedian = scoreFunctions[bottomMedianIndex].getScore();
            topWeight = sumWeights / 2 - cumulativeSumWeights[bottomMedianIndex];
        }

        return (topMedian * topWeight + bottomMedian * bottomWeight) / (topWeight + bottomWeight);
    }


    public double interdecileRange() {
        double out = interdecileRange;
        interdecileRange = -9999;
        return out;
    }

    private static class WeightComparator implements Comparator<NormalizedEnergyScore> {
        public int compare(NormalizedEnergyScore score1, NormalizedEnergyScore score2) {
            if (score1.getWeight() > score2.getWeight()) return 1;
            if (score1.getWeight() < score2.getWeight()) return -1;
            return 0;
        }
    }
}

