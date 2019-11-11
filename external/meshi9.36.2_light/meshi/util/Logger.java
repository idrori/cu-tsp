/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.*;
import meshi.optimizers.MCM;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;

import java.io.IOException;
import java.util.ArrayList;

public interface Logger {
    public abstract void mcm(ArrayList<Score> scoreFunctions, TotalEnergy energy,
                             int step,MCM.mcmStepResult mcmStepResult) throws IOException, UpdateableException, EvaluationException, AlignmentException;
    public double rms() ;
}