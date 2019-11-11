/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.util.*;

public class SummaAttribute implements MeshiAttribute {

    protected double fx, fy, fz;

    public int key() {
        return SUMMA_ATTRIBUTE;
    }

    public void reset() {
        fx = 0;
        fy = 0;
        fz = 0;
    }

}
