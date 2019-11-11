/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms;

import meshi.energy.*;

public abstract class NonBondedEnergyElement extends EnergyElement {
    protected double weight;
    public NonBondedEnergyElement() {
        super();
    }

    public abstract void set(Object obj);
    public void scaleWeight(double factor) {weight *= factor;}

}
