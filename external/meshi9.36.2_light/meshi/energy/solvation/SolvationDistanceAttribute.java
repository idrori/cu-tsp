/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

import meshi.geometry.Distance;
import meshi.util.*;
import meshi.util.mathTools.Sigma;

/**
 * In the SolvationEnergy evaluation function these fields are to be recalculated for EVERY
 * distance in the non-bonded list. They are used several times in the solvation evaluation
 * process, and we would like to calculate them once. To this end we attach this special
 * class as an attribute isOn each Distance instance.
 * For each distance we calculate two sigmoid values, and an hydrogen-bond strength value:
 * sigmCa1 - The carbon index of atom 2 isOn atom 1. If atom 2 is not a carbon then this value
 * should be zero. If atom 2 is a carbon then this value should be ~1.0 if atom 2
 * is spatially near atom 1. This index drops sigmoidally to zero the farther
 * atom 2 is.
 * sigmCa2 - The same as sigmCa1, except detailing the affect of atom 1 isOn atom 2.
 * sigmHBa1 - The hydrogen bond (HB) strength between atom 1 and 2. If atom 2 can not form HB
 * with atom 1, because of its chemical type then this value should be zero.
 * If atom 2 can form HB with atom 1 then this value should be ~1.0 if atom 1 and 2
 * are sufficiently close to create a HB, and if their orientation (defined
 * also by their base getAtoms - see HydrogenBondDahiyat) permit hydrogen bonding.
 * This value drops steeply to zero if the conditions to hydrogen bonding
 * are violated.
 * <p/>
 * Also provided are the carbon sigmoid values derivative relatives to the atom coordinates. They
 * have the general form: dsigmCa{1/2}d{x/y/z}{1/2}
 */


public class SolvationDistanceAttribute implements MeshiAttribute {


    public final Sigma sigma_12;
    public final Sigma sigma_21;
    public double sigmCa12;
    public double dsigmCa12dx;
    public double dsigmCa12dy;
    public double dsigmCa12dz;
    public double sigmCa21;
    public double dsigmCa21dx;
    public double dsigmCa21dy;
    public double dsigmCa21dz;

    public final Distance distance;
    public final SolvationDistanceType type;

    public SolvationDistanceAttribute(Distance distance, SolvationDistanceType type) {
        this(distance,type, null,null);
    }
    public SolvationDistanceAttribute(Distance distance, SolvationDistanceType type, SigmaParameters sigmaParameters_12, SigmaParameters sigmaParameters_21) {
        this.distance = distance;
        this.type = type;
        if (sigmaParameters_12 != null) {
            sigma_12 = new Sigma(sigmaParameters_12.end,
                    sigmaParameters_12.p1,
                    sigmaParameters_12.p2,
                    sigmaParameters_12.valAtp1,
                    sigmaParameters_12.valAtp2);
        } else sigma_12 = null;
        if (sigmaParameters_21 != null) {
            sigma_21 = new Sigma(sigmaParameters_21.end,
                    sigmaParameters_21.p1,
                    sigmaParameters_21.p2,
                    sigmaParameters_21.valAtp1,
                    sigmaParameters_21.valAtp2);
        } else sigma_21 = null;
    }


    public int key() {
        return SOLVATION_ALL_ATOM_ATTRIBUTE;
    }

    public void update() {
        double sigmas_tag;
        double dx = distance.dDistanceDx(), dy = distance.dDistanceDy(), dz = distance.dDistanceDz();
        double dis = distance.distance();
        // Handling Carbons
        // ----------------
        // Is atom2 non-polar? If so it should contribute to atom1's CNC
        sigma_12.sigma(dis);
        sigmCa12    = sigma_12.s();
        sigmas_tag  = sigma_12.s_tag();
        dsigmCa12dx = sigmas_tag * dx;
        dsigmCa12dy = sigmas_tag * dy;
        dsigmCa12dz = sigmas_tag * dz;

        // Is atom1 non-polar? If so it should contribute to atom2's CNC
        sigma_21.sigma(dis);
        sigmCa21    = sigma_21.s();
        sigmas_tag  = sigma_21.s_tag();
        dsigmCa21dx = sigmas_tag * dx;
        dsigmCa21dy = sigmas_tag * dy;
        dsigmCa21dz = sigmas_tag * dz;
    } // of update

}












