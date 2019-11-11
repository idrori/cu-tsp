/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

/**
 * All the atom type specific parameters of a solvation sigma function
 */
public class SigmaParameters {
    public final double end;
    public final double p1;
    public final double p2;
    public final double valAtp1;
    public final double valAtp2;


    public SigmaParameters(double end,
                           double p1,
                           double p2,
                           double valAtp1,
                           double valAtp2) {
        this.end = end;
        this.p1 = p1;
        this.p2 = p2;
        this.valAtp1 = valAtp1;
        this.valAtp2 = valAtp2;
    }
}
