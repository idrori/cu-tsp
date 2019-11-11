/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.TotalEnergy;

import java.util.*;

public class OptimizerException extends Exception {
    private final TotalEnergy energy;

    public OptimizerException(String msg) {
        super(msg + " ** No energy defined for this exception **");
        energy = null;
    }

    public OptimizerException(OptimizerException oe, TotalEnergy energy) {
        super(oe.getMessage() + "\n" + oe.getStackTrace().toString());
        this.energy = energy;
    }

    public TotalEnergy energy() {
        return energy;
    }
}


