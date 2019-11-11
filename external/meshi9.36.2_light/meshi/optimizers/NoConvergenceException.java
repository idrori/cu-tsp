/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import java.util.*;

public class NoConvergenceException extends Exception {
    public NoConvergenceException(int nSteps) {
        super("Minimization did not converge after " + nSteps + " steps");
    }
}

