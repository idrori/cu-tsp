/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.geometry.Distance;
import meshi.util.mathTools.Sigma;


public class SolvateDistanceAttributePolarNonPolar extends SolvateDistanceAttribute {
    public final Distance distance;

    public final SolvateDistanceType type() {
        return SolvateDistanceType.POLAR_HYOPHOBIC;
    }
    //public final int atom1number, atom2number;

    public final SigmaParameters sigmaParameters_12;


    // Carbon sigmoids and related derivatives
    public double sigmCa1;      /*1*/
    public double dsigmCa1dx;  /*2*/
    public double dsigmCa1dy;/*2*/
    public double dsigmCa1dz;/*2*/
    private final Sigma sigma = new Sigma();

    public SolvateDistanceAttributePolarNonPolar(Distance dis,
                                                 SigmaParameters sigmaParameters_12) {

        this.distance = dis;
        this.sigmaParameters_12 = sigmaParameters_12;
        resetSigmVals();
    }

    public final void update() {
        double sigmas_tag;

        // Handling Carbons
        // ----------------
        // Is atom2 non-polar? If so it should contribute to atom1's CNC
        sigma.sigma(distance.distance(), sigmaParameters_12.end,
                sigmaParameters_12.p1,
                sigmaParameters_12.p2,
                sigmaParameters_12.valAtp1,
                sigmaParameters_12.valAtp2);
        sigmCa1 = sigma.s();
        sigmas_tag = sigma.s_tag();
        dsigmCa1dx = sigmas_tag * distance.dDistanceDx();
        dsigmCa1dy = sigmas_tag * distance.dDistanceDy();
        dsigmCa1dz = sigmas_tag * distance.dDistanceDz();
    } // of update


    public final void resetSigmVals() {
        sigmCa1 = 0.0;
        dsigmCa1dx = 0.0;
        dsigmCa1dy = 0.0;
        dsigmCa1dz = 0.0;
    }
}