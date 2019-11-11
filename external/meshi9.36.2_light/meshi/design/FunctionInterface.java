/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.design;

/**
 * A representation of unitary mathematical function
 */
public interface FunctionInterface {
    /**
     * The function's value at a given position
     * @param x the independent  variable
     * @return   The function's value at x
     */
    public double valueAt(double x);
}
