/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.*;
import meshi.sequences.AlignmentException;
import meshi.util.*;

import java.io.IOException;

/**
 * Minimize energy according to a given setResidue of coordinates and an energy function
 */

public abstract class Optimizer {
    public final TotalEnergy energy;
    public final int maxSteps;
    public final int reportEvery;
    public static final Terminator optimizerTerminator = new Terminator();

    public enum OptimizerStatus {
        RUNNING, CONVERGED, UNCONVERGED, KILLED, DONE;
    }

    public Optimizer(TotalEnergy energy, int maxSteps, int reportEvery) {
        this.maxSteps = maxSteps;
        this.energy = energy;
        this.reportEvery = reportEvery;
        optimizerTerminator.reset();
    }

    public TotalEnergy energy() {
        return energy;
    }

    public abstract OptimizerStatus run() throws OptimizerException, UpdateableException, EvaluationException, AlignmentException, IOException;
}
