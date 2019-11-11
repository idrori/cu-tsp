/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.scoringFunctions;

import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.util.KeyWords;
import meshi.util.info.*;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 10/11/2010
 * Time: 23:18:00
 * To change this template use File | Settings | File Templates.
 */
public class NormalizedEnergyScore implements Comparable {
    private Sigmoid sigmoid;
    private Configuration configuration;
    private double score, weight;



    public NormalizedEnergyScore(Configuration configuration)  {
        sigmoid = new Sigmoid(configuration);
        this.configuration = configuration;

    }

    public void calcScore(ArrayList<MeshiInfo> energyInfo)  {
        double[] fieldValues      = getValues(energyInfo, configuration.fields);
        double[] weightedValues   = new double[configuration.fields.length];
        double[] normalizerValues = getValues(energyInfo, configuration.normalizers);
        double tempScore = 0;
        for (int iField = 0; iField < fieldValues.length; iField++) {
            double normalizationFactor = getNormalizationFactor(normalizerValues,configuration.exponentIndices[iField]);
            weightedValues[iField] = fieldValues[iField] *
                                     configuration.coefs[iField] *
                                     normalizationFactor;
            tempScore += weightedValues[iField];
            //System.out.println("yyyyyyy "+configuration.fields[iField]+" "+fieldValues[iField]+" "+configuration.coefs[iField]+" "+normalizationFactor+" "+tempScore);
        }
//        double tempScore = 0;
//        for (double d:weightedValues)
//            tempScore += d;
        score = sigmoid.backSigmoid(tempScore);
        weight = sigmoid.derivative(score);
    }

    public double getScore() {
        return score;
    }
    public double getWeight() {
        return weight;
    }


    private double getNormalizationFactor(double[] normalizersValues, int[] exponentsIndices){
        double factor = 1;
        for (int iNormalizer = 0; iNormalizer < normalizersValues.length; iNormalizer++){
            try {
                factor *= normalizersExponent(normalizersValues[iNormalizer],exponentsIndices[iNormalizer]);
            }
            catch (RuntimeException ex) {
                System.out.println("error in getNormalizationFactor\niNormalizer = "+iNormalizer);
                throw ex;
            }
        }
        return factor;
    }

    private double normalizersExponent(double normalizerValue, int exponentIndex){
        // Note that the indices here are taken from a MATLAB program and they (MATLAB people) never heard of index 0.
        if (normalizerValue < 0)
            throw new RuntimeException("Base cannot be smaller than 0, but it is "+normalizerValue);
        if (normalizerValue == 0) exponentIndex = 2;
        if (exponentIndex-1 >= configuration.exponentValues.length) return 0;
        if (configuration.exponentValues[exponentIndex-1] == 1)    return normalizerValue;
        if (configuration.exponentValues[exponentIndex-1] == 0.5)  return Math.sqrt(normalizerValue);
        if (configuration.exponentValues[exponentIndex-1] == 0)    return 1;
        if (configuration.exponentValues[exponentIndex-1] == -0.5) return 1/Math.sqrt(normalizerValue);
        if (configuration.exponentValues[exponentIndex-1] == -1)   return 1/normalizerValue;
        if (configuration.exponentValues[exponentIndex-1] == -2)   return 1/(normalizerValue*normalizerValue);
        if (configuration.exponentValues[exponentIndex-1] == -3)   return 1/(normalizerValue*normalizerValue*normalizerValue);
        throw new RuntimeException("weird value "+configuration.exponentValues[exponentIndex]);

    }

        private double[] getValues(ArrayList<MeshiInfo> energyInfo, String[] fields){
            double[] out = new double[fields.length];
            for (int iField = 0; iField < fields.length; iField++) {

                //*************** debug ***************
              // if (iField > 25) continue;

                //*************** debug ***************

                MeshiInfo element = null;
                for (MeshiInfo currentElement:energyInfo)  {
                    String tag = getTag(currentElement);
                    if (tag.equals(fields[iField])){
                        if (element != null)
                            throw new RuntimeException("field "+fields[iField]+" appears more than once "+element+" "+currentElement);
                        else
                            element = currentElement;
                    }
                }
                if (element == null)  {
//                    if (    (fields[iField].startsWith("cooperativeZstdS")) ||
//                            (fields[iField].startsWith("cooperativeStdS")))
//                        out[iField] = 0;
//                    else {
                        for (MeshiInfo currentElement : energyInfo)
                            System.out.println(getTag(currentElement));
                        throw new MissingFieldException("field " + fields[iField] + " is missing.");
 //                   }

                }
                else out[iField] = ((Double) element.getValue()).doubleValue();
            }

            return out;
        }
    private String getTag(MeshiInfo meshiInfo) {
        if ((meshiInfo.type == InfoType.WEIGHTED_MEDIAN_SCORE) ||
                (meshiInfo.type == InfoType.INTERDECILE))
            return meshiInfo.comment;
        return meshiInfo.type.tag;
    }

    public int compareTo(Object obj) {
        if (score - ((NormalizedEnergyScore) obj).score >0) return 1;
        if (score - ((NormalizedEnergyScore) obj).score <0) return -1;
        return 0;
    }
}
