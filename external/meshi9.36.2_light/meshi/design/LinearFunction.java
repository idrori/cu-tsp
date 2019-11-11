/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.design;

/**
 *  A representation of straight line.
 */
public class LinearFunction extends Function  {
    private double a,b;

    /**
     * @param a   slope
     * @param b   y-intercept
     */
    public LinearFunction(double a, double b) {
        this.a = a;
        this.b = b;
    }

    /**
     * The line's y-value for a given x-value
     * @param x independent variable
     * @return y - the function's value at x
     */
    public double valueAt(double x) {
        return a*x+b;
    }

    /**
     * The derivative of a straight line (its slope).
     * @return the slope.
     */
    protected FunctionInterface derived() {
        return new ConstantFunction(a);
    }
}
