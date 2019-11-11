/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.design;

/**
 * A  unitary mathematical function
 */
public abstract class Function implements FunctionInterface {

    /**
     *  @return The value of the function raised to the power of two.
     */
    public double square(double x) {
        double v = valueAt(x);
        return v*v;
    }

    /* *
    *    @return  The slope of the function at x
     */
    public double derivativeAt(double x){
        return derived().valueAt(x);
    }

    /**
     *
     * @return  The derivative of the function.
     */
    protected abstract FunctionInterface derived();
}
