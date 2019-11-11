 package meshi.design;
/**
 * A representation of sine wave.
 */
public class SinFunction extends TrigoFunction{

   /**
     * @param amplitude  is the peak deviation of the function from its center position.
     * @param frequency  specifies how many oscillations occur   per interval of the X-axis
     * @param phase specifies where in its cycle the oscillation begins
     */ public SinFunction(double amplitude, double frequency, double phase)  {
          super(amplitude, frequency, phase);
    }

    /**
    *    The basic cos function
     *  @param   x the x-axis
     *  @return cos(x)
     */
    protected  double  trigoElementaryFunction(double x) throws RuntimeException{
        System.out.println(Math.pow(-2,0.66666));
                return Math.sin(x);
    }

    /**
     *
     * @return return the derivative of sine, cosine.
     */
    protected FunctionInterface derived() {
        return new CosFunction(amplitude*frequency, frequency, phase);
    }

}
