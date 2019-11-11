/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.mathTools;

import java.util.*;

/**
 * The class creates a 1D cubic spline object from given parameters, and allows for the
 * calculation of the spline value and derivative at any given point inside the spline
 * intervals.
 * <p/>
 * After the method calc(x) is run, the public fields s and s_tag are updated:
 * s - the spline value at x.
 * s_tag - the spline first derivative at x.
 * <p/>
 * Note:
 * 1) calc will not work for values outside the range [break 1 , break n].
 * 2) The class operates isOn any setResidue of break points. However, if the break points are
 * evenly spaced then the calculation of calc(x) will be faster.
 * 3) See the constructor documention for more details isOn how to setResidue up the class properly.
 */

public class Spline1D {
    private double[] breaks;
    private double[][] coefs;
    private boolean evenBreaks = false;
    private double breakInterval;
    private int n;

    public double s = 0.0;
    public double s_tag = 0.0;

    /**
     * The constructor can create a specific spline object from a string with the format:
     * {break 1} {interval 1 X^3 coefficient} {interval 1 X^2 coef} {interval 1 X coef} {interval 1 constant coef} {break 2} {interval 2 X^3 coefficient} ... {break n}
     * <p/>
     * Note:
     * 1) There is an alternative constructor that requires a tokenizer of a string with
     * the above format.
     * 2) The breaks must increase monotonicly.
     * 3) The spline coefficients are not verified. Therefore, derivability is obtained only if the
     * coefficients are of a derivable spline.
     * 4) In order to have numerical stability, the relative precision of the values in the
     * dataLine must be at least 1e-10.
     */
    public Spline1D(String dataLine) {
        this(new StringTokenizer(dataLine));
    }

    public Spline1D(StringTokenizer st) {
        int i;
        boolean testEvenBreaks = true;
        if (((st.countTokens() - 1) / 5 != (st.countTokens() - 1) / 5.0) || (st.countTokens() < 6)) {
            throw new RuntimeException("Can not create a spline object from this line.Bad format");
        }
        breaks = new double[(st.countTokens() - 1) / 5 + 1];
        coefs = new double[(st.countTokens() - 1) / 5][4];
        for (i = 0; st.countTokens() > 1; i++) {
            breaks[i] = Double.valueOf(st.nextToken().trim()).doubleValue();
            coefs[i][0] = Double.valueOf(st.nextToken().trim()).doubleValue();
            coefs[i][1] = Double.valueOf(st.nextToken().trim()).doubleValue();
            coefs[i][2] = Double.valueOf(st.nextToken().trim()).doubleValue();
            coefs[i][3] = Double.valueOf(st.nextToken().trim()).doubleValue();
        }
        breaks[i] = Double.valueOf(st.nextToken().trim()).doubleValue(); // value of the last break point
        n = i;
        // testing for even breaks
        for (i = 1; i < breaks.length - 1; i++)
            if ((breaks[i + 1] - breaks[i]) != (breaks[1] - breaks[0]))
                testEvenBreaks = false;
        if (testEvenBreaks) {
            evenBreaks = true;
            breakInterval = breaks[1] - breaks[0];
        }
    }

    public final void calc(double x) {
        if ((x < breaks[0]) || (x > breaks[n]))
            throw new RuntimeException("X "+x+"is outside the spline range "+breaks[0]+"-"+breaks[n]);
        if (evenBreaks) {
            calcEvenBreaks(x);
            return;
        }
        int index;
        for (index = 0; (x < breaks[index]) || (x > breaks[index + 1]); index++) {
        }
        double offset = x - breaks[index];
        double offset2 = offset * offset;
        s = coefs[index][0] * offset2 * offset +
                coefs[index][1] * offset2 +
                coefs[index][2] * offset +
                coefs[index][3];

        s_tag = 3 * coefs[index][0] * offset2 +
                2 * coefs[index][1] * offset +
                coefs[index][2];
    }

    private final void calcEvenBreaks(double x) {
        int index = (int) ((x - breaks[0]) / breakInterval);
        double offset = x - breaks[0] - index * breakInterval;
        double offset2 = offset * offset;

        s = coefs[index][0] * offset2 * offset +
                coefs[index][1] * offset2 +
                coefs[index][2] * offset +
                coefs[index][3];

        s_tag = 3 * coefs[index][0] * offset2 +
                2 * coefs[index][1] * offset +
			coefs[index][2];		
	}

    public String toString() {
        String out =  "Spline1D evenBeaks = "+evenBreaks+"\tn = "+n+"\tbreakInterval = "+breakInterval+"\n";
        out += "s = "+s+"\ts_tag"+s_tag+"\n";
        out += "breaks:\n";
        for (double d: breaks) out += " "+d;
        out += "\ncoefs:\n";
        for (int i = 0; i < coefs.length; i++) {
            for (int j = 0; j < 4; j++)
                out += " "+coefs[i][j];
            out += "\n";
        }
        return out;
    }
}	
