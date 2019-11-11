/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.Parameters;
import meshi.parameters.AtomType;
import meshi.util.Utils;
import meshi.util.filters.Filter;

import java.util.StringTokenizer;

public class ExcludedVolParameters implements Parameters {
    public final double sigma;
    public final AtomType smallType;
    public final AtomType largeType;
    /**
     * width is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
     */
    public final double ALPHA1 = 1.0;
    public final double ALPHA2 = 0.2;
    /**
     * C is the parameter of  EV = C*(d-sigma)^4 in the range [0,sigma]
     */
    public final double C1, C2;


    public ExcludedVolParameters() {
        smallType = AtomType.XXX;
        largeType = AtomType.XXX;
        sigma = -1;
        C1 = -1;
        C2 = -1;
    }

    public ExcludedVolParameters(StringTokenizer st) {
        AtomType first = AtomType.type(st.nextToken());
        AtomType second = AtomType.type(st.nextToken());
        if (first.ordinal() > second.ordinal()) {
            smallType = second;
            largeType = first;
        } else {
            smallType = first;
            largeType = second;
        }
        sigma = Utils.toDouble(st.nextToken());
        C1 = 1 / (ALPHA1 * ALPHA1);
        C2 = 1 / (ALPHA2 * ALPHA2 * ALPHA2 * ALPHA2);
    }

    public Filter isA() {
        return (new isA());
    }

    public Parameters create(StringTokenizer stringTokenizer) {
        return (new ExcludedVolParameters(stringTokenizer));
    }

    private class isA implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof ExcludedVolParameters);
        }
    }

    public String toString() {
        return "" + sigma;
    }
}
