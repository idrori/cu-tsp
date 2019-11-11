/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.design;

/**
 *  A representation of a horizontal line.
 */
public class ConstantFunction extends LinearFunction {
    private double c;

    /**
     * @param c  y-intercept
     */
    public ConstantFunction(double c) {
        super(0,c);
     }

}
