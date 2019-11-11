/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.bond.BondEnergy;
import meshi.energy.simpleEnergyTerms.bond.BondEnergyElement;
import meshi.molecularElements.atoms.Atom;
import meshi.util.UpdateableException;

/**
 * This class provide a simple way to find the step length. It start with a certain length and check if
 * it brings to energy reduction. If it does not it is shortened by a reduction factor. This repeats until a
 * step size is found that produce energy reduction, or the step size is so small that an exception is thrown.
 * After the step size is returned, it is multiplied by an expansion factor, so that the next call to this
 * class evaluation is started with a longer guess. If the expension factor is setResidue to 0, the step size initial
 * guess is always 1 (as should be tried for all quasi-newton techniques).
 * <p/>
 * How to use this class:
 * ----------------------
 * a) Instantiate this class with the desired parameters.
 * b) In a position Xk: evaluate the energy gradients and coordinates to position Xk.
 * c) Run 'findStepLength(Vec[n][2])' where the first column in Vec is the position Xk, and the second
 * column is the direction Pk. This method returns the found step length. This method also changes the coordinates
 * in class 'energy' to the coordinates in the new (minimized) position: Xk+1 = Xk + (step_length)*Pk. Also the
 * gradients in the energy class are updated.
 * d) Check for thrown exceptions to make sure that the step length is correct.
 */

public class SimpleStepLength extends LineSearch {

    private double stepSize, stepSizeReduction, stepSizeExpansion;
    private double old_e, new_e;
    private double[][] coordinates;
    private final double verySmall = Math.exp(-60);
    private final double INFINITY = 1 / 0.0;
    int i;


    public SimpleStepLength(TotalEnergy energy,
                            double stepSize1,
                            double stepSizeReduction,
                            double stepSizeExpansion) {
        super(energy);
        this.stepSizeExpansion = stepSizeExpansion;
        this.stepSizeReduction = stepSizeReduction;
        coordinates = energy.coordinates();
        if (stepSizeExpansion > 0)
            this.stepSize = stepSize1 / stepSizeExpansion;
        else
            this.stepSize = 1;
        if (stepSizeReduction <= 0)
            throw new RuntimeException("\n\nIrrecoverable error in the line search method: SimpleStepLength.\n" +
                    "The step size reduction parameter in the constructor is non-positive.\n");
    }

    public double findStepLength(double[][] inputCoordinates) throws LineSearchException, UpdateableException, EvaluationException{
        double step;
        if (inputCoordinates == coordinates)
            throw new LineSearchException(LineSearchException.WEIRD_INPUT_TO_FIND_STEP_LENGTH,
                    "\n\nThe input array to the function 'findStepLength' " +
                            "has the same pointer as the 'coordinate' array in energy. \n" +
                            "It should be a different array.\n");
        stepSize = stepSize * stepSizeExpansion / stepSizeReduction;
        if (stepSizeExpansion <= 0)
            stepSize = 1;
        for (i = 0; i < coordinates.length; i++) {
            coordinates[i][0] = inputCoordinates[i][0];
        }
        old_e = energy.getLastEnergy();
        new_e = old_e;
        // If no energy reduction is achieved for a specific step size it is reduced
        // until a step that produce reduction in energy is found.
        while (new_e >= old_e) {
            stepSize *= stepSizeReduction;
            if (stepSize < verySmall) {
                throw new LineSearchException(LineSearchException.NOT_A_DESCENT_DIRECTION,
                        "\n\nStepSize = " + stepSize + "\nverySmall=" + verySmall + "\nstepSizeExpansion=" + stepSizeExpansion + "\n" +
                                "stepSizeReduction=" + stepSizeReduction + "\n" +
                                "new_e=" + new_e + "\n" +
                                "old_e=" + old_e + "\n" +
                                energy.report(-11111) +
                                "The search direction is apparently not a descent direction. \n" +
                                "This problem might be caused by incorrect diffrentiation " +
                                "of the energy function,\n" +
                                "or by numerical instabilities of the minimizing techniques " +
                                "(such as not fullfilling the Wolf condtions in BFGS).\n");
            }
            for (i = 0; i < coordinates.length; i++)     {
                step =  stepSize * inputCoordinates[i][1];
               if (step> 1) step = 1 ;
                if (step < -1) step = -1;
                coordinates[i][0] = inputCoordinates[i][0] + step;
                //if (coordinates[i][0] > 10000) throw new RuntimeException("This is weird "+i +"   "+coordinates[i][0] +"   "+step+"     "+ stepSize+"    "+inputCoordinates[i][1]);
            }
                new_e = energy.evaluate();   // The energy at the new coordinates.
        }
        return stepSize;
    }


}    
