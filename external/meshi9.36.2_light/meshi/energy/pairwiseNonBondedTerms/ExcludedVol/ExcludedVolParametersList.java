/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.Parameters;
import meshi.geometry.Distance;
import meshi.parameters.AtomType;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;

import java.util.ArrayList;
import java.util.StringTokenizer;


public class ExcludedVolParametersList extends ArrayList<ExcludedVolParameters> {

    public ExcludedVolParametersList(String parametersFileName) {
        super();

        ExcludedVolParameters tmpParam;
        AtomType maxType = AtomType.XXX;
        ExcludedVolParameters[] tmpAr;
        MeshiLineReader lines;
        String line;

        if (parametersFileName == null)
            throw new RuntimeException("No parameter file name in " + this);

        Utils.println("Loading " + this + " parameters from " + parametersFileName);
        // First reading to see that everything is ok, and also to see what is the
        // largest atom type.
        try {
            lines = new MeshiLineReader(parametersFileName);
            while ((line = lines.readLine("#")) != null) {
                tmpParam = (ExcludedVolParameters) createParameters(line);
                if (tmpParam.largeType.ordinal() > maxType.ordinal())
                    maxType = tmpParam.largeType;
            }
        }
        catch (RuntimeException e) {
            System.out.println("A problem while reading parameters file " +
                    parametersFileName);
            throw e;
        }
        tmpAr = new ExcludedVolParameters[((maxType.ordinal() + 1) * (maxType.ordinal() + 2)) / 2];
        // Second reading actually creates the ExcludedVolParameters objects.
        try {
            lines = new MeshiLineReader(parametersFileName);
            while ((line = lines.readLine("#")) != null) {
                tmpParam = (ExcludedVolParameters) createParameters(line);
                tmpAr[tmpParam.largeType.ordinal() * (tmpParam.largeType.ordinal() + 1) / 2 + tmpParam.smallType.ordinal()] = tmpParam;
            }
        }
        catch (RuntimeException e) {
            System.out.println("A problem while reading parameters file " +
                    parametersFileName);
            throw e;
        }

        for (ExcludedVolParameters aTmpAr : tmpAr) add(aTmpAr);
    }

    public Parameters createParameters(String line) {
        return new ExcludedVolParameters(new StringTokenizer(line));
    }

    public Parameters parameters(Object obj) {
        Distance distance = (Distance) obj;
        AtomType largeType = distance.largeType;
        AtomType smallType = distance.smallType;
        try {
            return get(largeType.ordinal() * (largeType.ordinal() + 1) / 2 + smallType.ordinal());
        }
        catch (Exception e) {
            throw new RuntimeException("ExcludedVolParametersList: " + largeType + " " +
                    smallType + "\n" + e);
        }
    }
}
