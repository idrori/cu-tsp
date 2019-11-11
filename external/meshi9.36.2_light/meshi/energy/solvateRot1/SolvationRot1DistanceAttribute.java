/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvateRot1;

import meshi.energy.solvation.SigmaParameters;
import meshi.energy.solvation.SolvationDistanceAttribute;
import meshi.geometry.Distance;

/**
 * In the SolvateRot1Energy evaluation function these fields are to be recalculated for EVERY
 * distance in the non-bonded list. They are used several times in the solvation evaluation
 * process, and we would like to calculate them once. To this end we attach this special
 * class as an attribute isOn each Distance instance.
 * For each distance we calculate two sigmoid values: one for atom 1 in the distance (a1),
 * and one for atom 2 in the distance (a2):
 * sigmCa1 - The carbon index of atom 2 isOn atom 1. This value should be ~1.0 if atom 2
 * is spatially near atom 1. This index drops sigmoidally to zero the farther
 * atom 2 is.
 * sigmCa2 - The same as sigmCa1, except detailing the affect of atom 1 isOn atom 2.
 * <p/>
 * Also provided are the sigmoid values derivative relatives to the atom coordinates. They
 * have the general form: dsigmCa{1/2}d{x/y/z}{1/2}
 */


public class SolvationRot1DistanceAttribute extends SolvationDistanceAttribute {

    public SolvationRot1DistanceAttribute(Distance distance1, SigmaParameters sigmaParameters_12, SigmaParameters sigmaParameters_21) {
        super(distance1, null, sigmaParameters_12,sigmaParameters_21);
    }

    public int key() {
        return SOLVATE_ROT1_ATTRIBUTE;
    }

    public double sigmCa1;
    public double sigmCa2;
    public double dsigmCa1dx1;
    public double dsigmCa1dy1;
    public double dsigmCa1dz1;
    public double dsigmCa1dx2;
    public double dsigmCa1dy2;
    public double dsigmCa1dz2;
    public double dsigmCa2dx1;
    public double dsigmCa2dy1;
    public double dsigmCa2dz1;
    public double dsigmCa2dx2;
    public double dsigmCa2dy2;
    public double dsigmCa2dz2;


    public final void resetAllSigmVals() {
        sigmCa1 = 0.0;
        dsigmCa1dx1 = 0.0;
        dsigmCa1dy1 = 0.0;
        dsigmCa1dz1 = 0.0;
        dsigmCa1dx2 = 0.0;
        dsigmCa1dy2 = 0.0;
        dsigmCa1dz2 = 0.0;
        sigmCa2 = 0.0;
        dsigmCa2dx1 = 0.0;
        dsigmCa2dy1 = 0.0;
        dsigmCa2dz1 = 0.0;
        dsigmCa2dx2 = 0.0;
        dsigmCa2dy2 = 0.0;
        dsigmCa2dz2 = 0.0;
    }

}