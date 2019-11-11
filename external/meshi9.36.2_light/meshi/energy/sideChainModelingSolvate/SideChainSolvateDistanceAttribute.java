/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.sideChainModelingSolvate;

import meshi.util.*;

/**
 * This class is a truncated form of the corresponding class in the "meshi.energy.solvation" package. It was designed
 * solely for accelerating the application SCMOD (concurrent sidechain modeling), and should not be used for other
 * purposes.
 * <p/>
 * Please do not use it!
 */


public class SideChainSolvateDistanceAttribute implements MeshiAttribute {

    public SideChainSolvateDistanceAttribute() {
    }

    public int key() {
        return SIDE_CHAIN_SOLVATE_ALL_ATOM_ATTRIBUTE;
    }

    public double sigmCa1;
    public double sigmHBa1;
    public double sigmCa2;
    public double sigmHBa2;

    private final void resetHBSigmVals() {
        sigmHBa1 = 0.0;
        sigmHBa2 = 0.0;
    }


    public final void resetAllSigmVals() {
        sigmCa1 = 0.0;
        sigmCa2 = 0.0;
        resetHBSigmVals();
    }


}
