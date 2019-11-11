package meshi.design;
/**
 * A master class for all trigonometric functions
 */
public abstract class TrigoFunction extends Function {
    protected double amplitude, frequency, phase;

    /**
     * @param amplitude  is the peak deviation of the function from its center position.
     * @param frequency  specifies how many oscillations occure   per interval of the X-axis
     * @param phase pecifies where in its cycle the oscillation begins
     */
    public TrigoFunction(double amplitude, double frequency, double phase) {
        this.amplitude = amplitude;
        this.frequency = frequency;
        this.phase = phase;
    }

    /**
     * The value of the trigonometric function at a given position.
     * @param x independent variable
     * @return    the function's value at x.
     */
    public double valueAt(double x) {
           return amplitude*trigoElementaryFunction(frequency*x+phase);
       }

    /**
     * Concrete trigonometric function
     * @param x    independent variable
     * @return      the function's value at x.
     */
       protected abstract double  trigoElementaryFunction(double x);
    
}
