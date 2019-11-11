/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.*;
import meshi.util.UpdateableException;

/**
 * This class implements a line search that satisfies the Wolf conditions according to the scheme in: Numerical Optimization
 * by J. Nocendal & S. J. Wright, Springer 1999, pp 55-62.
 * <p/>
 * Being in position Xk, and having a descent search direction Pk, a step length (width) is sought, that leads to a new position Xk+1,
 * so that Xk+1 = Xk + width*Pk, and stisfies the 2 Wolf conditons:
 * 1) E(Xk+1) < E(Xk) + c1*width*grad(Xk)*Pk   (sufficient decrease)
 * 2) |grad(Xk+1)*Pk| < c2*|grad(Xk)*Pk|   (curvature condition)
 * The implementation follows the one in the above reference. In short, the variable width is iteratively enlarged, until it
 * brackets (with the previous iteration value) an interval that contain some points with the Wolf conditions. The function
 * Zoom is then called, that finds the exact width where the Wolf conditions hold. The initial choice of width (at the first
 * iteration) is the popular min{1 , -2*(E(Xk-1) - E(Xk))/(grad(Xk)*Pk)}. Other choices are possible to code if this one
 * leads to failure of the line search algorithm. This algorithm is safe guarded against infinite loop, by a condition isOn the
 * maximal number of energy evaluations. See bellow what to do if this condition is violated.
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
 * <p/>
 * <p/>
 * <p/>
 * Parameters (Good default values are given after the parameter name)
 * ----------
 * c1,c2 - 1e-4,0.9 - The two paramters of the Wolf conditions. Must satisfy: 0<c1<c2<1
 * maxNumEvaluations - 10 - The maximal number of step length trails allowed. This gives an upper limit isOn the total number of
 * evaluations both in the bracketing and Zoom. If line search fails because this number was
 * exceeded try to enlarge 'extendAlphaFactor' or change the initial width guess. This error might
 * also be caused if the tolerance of the minimizer is setResidue to be extremely small.
 * extendAlphaFactor - 3 - After a certain step length fails, the next step length is taken as this multiple of the failing
 * last width.
 * <p/>
 * Constant Parameters
 * -------------------
 * interSafeGuardFactor - 20 - This factor is used in the calculation of the cubic interpolation safeguard. It's exact value
 * is not crucial to the minimizer performance. 20 is generally good.
 */

public class WolfConditionLineSearch extends LineSearch {

    private double maxNumEvaluations;
    private double extendAlphaFactor;
    private double c1, c2;
    private final double interSafeGuardFactor = 20;
    private double interSafeGuard;
    private double oldE, newE;
    private double[][] coordinates;
    private int firstRun;
    private double alphaI, alphaI1, alphaFinal;
    private int n; // number of variables
    private double e0, eI, eI1; //Energy values along Pk
    private double grad0, gradI, gradI1; //gradients along Pk
    private int numAlphaEvaluations, stop;
    private int i;


    public WolfConditionLineSearch(TotalEnergy energy,
                                   double c1, double c2,
                                   double extendAlphaFactor,
                                   int maxNumEvaluations) {
        super(energy);
        this.c1 = c1;
        this.c2 = c2;
        this.maxNumEvaluations = maxNumEvaluations;
        this.extendAlphaFactor = extendAlphaFactor;
        interSafeGuard = extendAlphaFactor / interSafeGuardFactor;
        coordinates = energy.coordinates();
        n = coordinates.length;
        Reset();
    }

    public void Reset() {
        firstRun = 1;
        alphaI1 = 1.0;
    }

    public void Reset(double resAlpha) {
        firstRun = 1;
        alphaI1 = resAlpha;
    }

    public double findStepLength(double[][] inputCoordinates) throws LineSearchException, UpdateableException, EvaluationException {
        if (inputCoordinates == coordinates)
            throw new LineSearchException(LineSearchException.WEIRD_INPUT_TO_FIND_STEP_LENGTH,
                    "\n\nThe input array to the function 'findStepLength' has the same pointer" +
                            " as the 'coordinate' array in energy. \n" +
                            "It should be a different array.\n");
        //Initializing the run
        newE = energy.getLastEnergy();
        grad0 = 0; // calculating the gradient at a=0
        for (i = 0; i < n; i++)
            grad0 -= inputCoordinates[i][1] * coordinates[i][1];
        if (grad0 >= 0)
            throw new LineSearchException(LineSearchException.NOT_A_DESCENT_DIRECTION,
                    "\n\nThe search direction is not a descent direction. \n" +
                            "This problem might be caused by incorrect diffrentiation of the energy function,\n" +
                            "or by numerical instabilities of the minimizing techniques" +
                            "(such as not fullfilling the Wolf condtions in BFGS).\n");
        numAlphaEvaluations = -1;
        stop = 0;
        e0 = newE;
        eI = e0;
        alphaI = 0;
        // Choosing the initial width guess
        if (firstRun != 0) {
            firstRun = 0;
        } else {
            alphaI1 = 2.02 * (newE - oldE) / grad0;
            if (alphaI1 > 1)
                alphaI1 = 1;
        }
        oldE = newE; // For the next time line search is called
        // Bracketing the Wolf area
        while ((numAlphaEvaluations <= maxNumEvaluations) && (stop == 0)) {
            numAlphaEvaluations++;
            for (i = 0; i < n; i++)
                coordinates[i][0] = inputCoordinates[i][0] + alphaI1 * inputCoordinates[i][1];
            eI1 = energy.evaluate();
            gradI1 = 0; // calculating the gradient at a=alphaI1
            for (i = 0; i < n; i++)
                gradI1 -= inputCoordinates[i][1] * coordinates[i][1];
            if ((eI1 > (e0 + c1 * alphaI1 * grad0)) || ((eI1 >= eI) && (numAlphaEvaluations > 0))) {
                Zoom(inputCoordinates, 0); //calling the regular zoom
            } else {
                if (Math.abs(gradI1) <= (-c2 * grad0)) {
                    alphaFinal = alphaI1;
                    stop = 1;
                } else {
                    if (gradI1 >= 0) {
                        Zoom(inputCoordinates, 1); //calling the inverse Zoom
                    }
                }
            }
            alphaI = alphaI1;
            gradI = gradI1;
            eI = eI1;
            alphaI1 = alphaI1 * extendAlphaFactor;
        }
        if (numAlphaEvaluations > maxNumEvaluations) {
            for (i = 0; i < n; i++)  // Returining the coordinates to the original state
                coordinates[i][0] = inputCoordinates[i][0];
            energy.evaluate();
            throw new LineSearchException(LineSearchException.WOLF_CONDITION_NOT_MET,
                    "\n\nWolf conditions not met. The line search did not converge,\n" +
                            "and exceeded the maximal number of step extensions allowed. See the help\n" +
                            "isOn this class for possible solutions to this problem.\n");
        }
        return alphaFinal;
    }

    // The function Zoom finds a step length satisfing the Wolf conditions, given the bracketing of alphaI and alphaI1.
    // It was separated into a different function to make the code more readable.

    private void Zoom(double[][] inputCoordinates, int inv) throws UpdateableException, EvaluationException {
        double alphaHi, alphaLow, alphaNew = 0;
        double eHi, eLow, eNew = 0;
        double gradHi, gradLow, gradNew = 0;
        double a = 0, b = 0, ga = 0, gb = 0, ea, eb, interval, d1 = 0, d2 = 0; // Cubic interpolation variables

        eHi = eI1;
        gradHi = gradI1;
        alphaHi = alphaI1;
        eLow = eI;
        gradLow = gradI;
        alphaLow = alphaI;
        if (inv == 1) {
            eHi = eI;
            gradHi = gradI;
            alphaHi = alphaI;
            eLow = eI1;
            gradLow = gradI1;
            alphaLow = alphaI1;
        }

        while ((numAlphaEvaluations <= maxNumEvaluations) && (stop == 0)) {
            numAlphaEvaluations++;
            // Cubic interpolation of the next step length
            a = alphaLow;
            b = alphaHi;
            ga = gradLow;
            gb = gradHi;
            ea = eLow;
            eb = eHi;
            if (a > b) { //switch
                b = alphaLow;
                a = alphaHi;
                gb = gradLow;
                ga = gradHi;
                eb = eLow;
                ea = eHi;
            }
            d1 = ga + gb - 3 * (ea - eb) / (a - b);
            if ((d1 * d1 - ga * gb) >= 0) {
                d2 = Math.sqrt(d1 * d1 - ga * gb);
                alphaNew = b - (b - a) * (gb + d2 - d1) / (gb - ga + 2 * d2);
            } else
                alphaNew = a; // Forcing bisection
            interval = Math.abs(b - a);
            if ((Math.abs(a - alphaNew) < (interval * interSafeGuard)) || (Math.abs(b - alphaNew) < (interval * interSafeGuard)))
                alphaNew = (a + b) / 2; //If the cubic minimum is to close to the edges, the bisection method is choosen
            // Continue with zoom - calculating the properties of the new found step length
            for (i = 0; i < n; i++)
                coordinates[i][0] = inputCoordinates[i][0] + alphaNew * inputCoordinates[i][1];
            eNew = energy.evaluate();
            gradNew = 0; // calculating the gradient at a=alphaNew
            for (i = 0; i < n; i++)
                gradNew -= inputCoordinates[i][1] * coordinates[i][1];
            if ((eNew > (e0 + c1 * alphaNew * grad0)) || (eNew >= eLow)) {
                alphaHi = alphaNew;
                eHi = eNew;
                gradHi = gradNew;
            } else {
                if (Math.abs(gradNew) <= (-c2 * grad0))
                    stop = 1;
                else {
                    if ((gradNew * (alphaHi - alphaLow)) >= 0) {
                        alphaHi = alphaLow;
                        eHi = eLow;
                        gradHi = gradLow;
                    }
                    alphaLow = alphaNew;
                    eLow = eNew;
                    gradLow = gradNew;
                }
            }
        }
        alphaFinal = alphaNew;
    }

}    
