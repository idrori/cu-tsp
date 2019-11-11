/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.design;

/**
 * A representation of cosine wave.
 */
public class CosFunction extends TrigoFunction{
    /**
     * @param amplitude  is the peak deviation of the function from its center position.
     * @param frequency  specifies how many oscillations occur   per interval of the X-axis
     * @param phase specifies where in its cycle the oscillation begins
     */
    public CosFunction(double amplitude, double frequency, double phase)  {
          super(amplitude, frequency, phase);
    }

    /**
    *    The basic cos function
     *  @param   x the x-axis
     *  @return cos(x)
     */
    public   double  trigoElementaryFunction(double x) {
            return Math.cos(x);
    }

    /**
     *
     * @return return the derivative of cosine, -sin.
     */
    protected FunctionInterface derived() {
        return new SinFunction(-amplitude*frequency, frequency, phase);
    }
}
