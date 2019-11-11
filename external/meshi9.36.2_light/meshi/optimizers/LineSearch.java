/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.*;
import meshi.util.UpdateableException;

/**
 * An abstract class of the common function calls in all the Line-Search procedure.
 * <p/>
 * If the simulation is at point Xk,and a descent direction is Pk, findStepLength gives a scalar Ak,
 * that says how much to travel along Pk. So that Xk+1 = Xk + Ak*Pk. The coordinates of the initial point (Xk)
 * should be given in the first column of the argument 'coordinates'. The direction of search (Pk)
 * (direction + magnitude) from the initial point sould be given in the second column of 'coordinates'.
 */

public abstract class LineSearch {

    protected TotalEnergy energy;

    public LineSearch(TotalEnergy energy) {
        this.energy = energy;
    }

    public abstract double findStepLength(double[][] coordinates) throws LineSearchException, UpdateableException, EvaluationException;
}
