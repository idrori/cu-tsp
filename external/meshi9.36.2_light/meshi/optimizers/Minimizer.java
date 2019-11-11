/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.sequences.AlignmentException;
import meshi.util.Terminator;
import meshi.util.UpdateableException;
import meshi.util.Utils;

import java.io.IOException;

/**
 * Minimize energy according to a given setResidue of coordinates and an energy function
 */

public abstract class Minimizer extends Optimizer {
    public final int MAX_KICKSTARTS = 10;
    public double tolerance;
    private double forceMagnitude;
    private int numberOfKickStrarts;
    public static final Terminator terminator = new Terminator();

    public Minimizer(TotalEnergy energy, int maxSteps, int reportEvery, double tolerance) throws UpdateableException {
        super(energy, maxSteps, reportEvery);
        this.tolerance = tolerance;
    }

    public OptimizerStatus run() throws OptimizerException, UpdateableException, EvaluationException, AlignmentException, IOException {
        return run(true);
    }

    public OptimizerStatus run(boolean testFlag) throws OptimizerException, UpdateableException, EvaluationException {
        OptimizerStatus os = init();
        if (os == OptimizerStatus.CONVERGED) {
            Utils.println("*************** Gradient is too low to minimize *******");
            double e = energy.evaluate();
            double grad = energy().getGradMagnitude();
            Utils.println("Energy is "+e);
            Utils.println("Gradient  is "+grad);
            return OptimizerStatus.CONVERGED;
        }
        numberOfKickStrarts = 0;
        int step;
        boolean minimizationStepOK;
        for (step = 1; status(step) == OptimizerStatus.RUNNING; step++) {
            try {
                minimizationStepOK = minimizationStep();
            } catch (Exception e) {
                System.out.println("Error in minimization step "+step);
                e.printStackTrace();
                System.out.println("---------------------------------------------------------------------------------------------------");
                throw new RuntimeException(e);
            }
            if (!minimizationStepOK) {
                if (numberOfKickStrarts >= MAX_KICKSTARTS)
                    throw new OptimizerException("\n\nThe simulation was restarted for " + MAX_KICKSTARTS + " times " +
                            "which is more than allowed.\n" +
                            "So many restarts are indicative of an ill-shaped energy function or " +
                            "an energy differentiation\n");
                try {
                    kickStart();
                    Utils.println("kickstart # " + numberOfKickStrarts + " done");
                }
                catch (OptimizerException oe) {
                    if ((testFlag) && Utils.verbose()) energy.test();
                    throw oe;
                }
                numberOfKickStrarts++;
            }

            if (step % reportEvery == 0)
                Utils.println(energy().report(step));
        }
        return status(step);
    }

    public OptimizerStatus status(int step) {
        if (terminator.dead()) {
            return OptimizerStatus.KILLED;
        }
        forceMagnitude = energy.getGradMagnitude();

        if (forceMagnitude < tolerance) return OptimizerStatus.CONVERGED;
        if (step <= maxSteps) return OptimizerStatus.RUNNING;
        return OptimizerStatus.UNCONVERGED;
    }

    protected abstract OptimizerStatus init() throws OptimizerException, UpdateableException, EvaluationException;

    protected abstract boolean minimizationStep() throws OptimizerException, UpdateableException, EvaluationException;

    protected abstract OptimizerStatus kickStart() throws OptimizerException, UpdateableException, EvaluationException;
}
