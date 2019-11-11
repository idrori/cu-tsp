/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.mathTools;

import meshi.util.*;
import meshi.util.file.MeshiWriter;

import java.util.*;

/**
 * The class creates a 2D cubic spline object from given parameters, and allows for the
 * calculation of the spline value and derivative at any given point inside the spline
 * intervals.
 * <p/>
 * After the method calc(x,y) is run, the public fields s s_tag_x and s_tag_y are updated:
 * s - the spline value at x,y.
 * s_tag_x - the spline first derivative at x,y with respect to x.
 * s_tag_y - the spline first derivative at x,y with respect to y.
 * <p/>
 * Note:
 * 1) calc will not work for values outside the range [break 1 , break n].
 * 2) The class operates isOn any setResidue of break points. However, if the break points are
 * evenly spaced (in both axes) then the calculation of calc(x,y) will be faster.
 * 3) See the constructor documention for more details isOn how to setResidue up the class properly.
 */

public class Spline2D {
    private double[] breaksX;
    private double[] breaksY;
    private double[][][] coefs;
    private boolean evenBreaks = false;
    private double breakIntervalX;
    private double breakIntervalY;

    public double s = 0.0;
    public double s_tag_x = 0.0;
    public double s_tag_y = 0.0;
    private static double testDiff = 0;

    /**
     * The constructor can create a specific spline object from a string with the format (Note that following
     * 6 lines should appear in the SAME line, but I broke it for clarity). Just to make it perfectly clearRows,
     * this object is created from ONE very long test-line:
     * {Number of breaks in X axis} {Number of breaks in Y axis} {break 1X} ... {break nX} {break 1Y} ... {break nY}
     * {square(1,1) constant coef} {square(1,1) X coef} {square(1,1) X^2 coef} {square(1,1) X^3 coefficient}
     * {square(1,1) Y coef} {square(1,1) Y*X coef} {square(1,1) Y*X^2 coef} {square(1,1) Y*X^3 coefficient}
     * {square(1,1) Y^2 coef} {square(1,1) Y^2*X coef} {square(1,1) Y^2*X^2 coef} {square(1,1) Y^2*X^3 coefficient}
     * {square(1,1) Y^3 coef} {square(1,1) Y^3*X coef} {square(1,1) Y^3*X^2 coef} {square(1,1) Y^3*X^3 coefficient}
     * {square(1,2) constant coef} {square(1,2) X coef} ...
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
    public Spline2D(String dataLine) {
        this(new StringTokenizer(dataLine));
    }

    public Spline2D(StringTokenizer st) {
        int i , j, k;
        int numXbreak = Integer.valueOf(st.nextToken().trim()).intValue();
        int numYbreak = Integer.valueOf(st.nextToken().trim()).intValue();
        breaksX = new double[numXbreak];
        breaksY = new double[numYbreak];
        coefs = new double[numXbreak][numYbreak][16];
        for (i = 0; i < numXbreak; i++)
            breaksX[i] = Double.valueOf(st.nextToken().trim()).doubleValue();
        for (i = 0; i < numYbreak; i++)
            breaksY[i] = Double.valueOf(st.nextToken().trim()).doubleValue();
        for (i = 0; i < (numXbreak - 1); i++)
            for (j = 0; j < (numYbreak - 1); j++)
                for (k = 0; k < 16; k++)
                    coefs[i][j][k] = Double.valueOf(st.nextToken().trim()).doubleValue();
        // testing for even breaks
        evenBreaks = testForEvenBreaks(breaksX) && testForEvenBreaks(breaksY);
        if (evenBreaks) {
            breakIntervalX = breaksX[1] - breaksX[0];
            breakIntervalY = breaksY[1] - breaksY[0];
        }
    }


    public final boolean calc(double x, double y) {
        if (evenBreaks) {
            return calcEvenBreaks(x, y);
        }
        int indexX;
        int indexY;
        for (indexX = 0; (x > breaksX[indexX + 1]); indexX++) {
        }
        for (indexY = 0; (y > breaksY[indexY + 1]); indexY++) {
        }
        if ((indexX >= coefs.length) || (indexY >= coefs[0].length)) return false;
        double offsetX = x - breaksX[indexX];
        double offsetX2 = offsetX * offsetX;
        double offsetX3 = offsetX2 * offsetX;
        double offsetY = y - breaksY[indexY];
        double offsetY2 = offsetY * offsetY;
        double offsetY3 = offsetY2 * offsetY;
        s = coefs[indexX][indexY][0] +
                coefs[indexX][indexY][1] * offsetX +
                coefs[indexX][indexY][2] * offsetX2 +
                coefs[indexX][indexY][3] * offsetX3 +
                coefs[indexX][indexY][4] * offsetY +
                coefs[indexX][indexY][5] * offsetY * offsetX +
                coefs[indexX][indexY][6] * offsetY * offsetX2 +
                coefs[indexX][indexY][7] * offsetY * offsetX3 +
                coefs[indexX][indexY][8] * offsetY2 +
                coefs[indexX][indexY][9] * offsetY2 * offsetX +
                coefs[indexX][indexY][10] * offsetY2 * offsetX2 +
                coefs[indexX][indexY][11] * offsetY2 * offsetX3 +
                coefs[indexX][indexY][12] * offsetY3 +
                coefs[indexX][indexY][13] * offsetY3 * offsetX +
                coefs[indexX][indexY][14] * offsetY3 * offsetX2 +
                coefs[indexX][indexY][15] * offsetY3 * offsetX3;

        s_tag_x = coefs[indexX][indexY][1] +
                coefs[indexX][indexY][2] * 2 * offsetX +
                coefs[indexX][indexY][3] * 3 * offsetX2 +
                coefs[indexX][indexY][5] * offsetY +
                coefs[indexX][indexY][6] * offsetY * 2 * offsetX +
                coefs[indexX][indexY][7] * offsetY * 3 * offsetX2 +
                coefs[indexX][indexY][9] * offsetY2 +
                coefs[indexX][indexY][10] * offsetY2 * 2 * offsetX +
                coefs[indexX][indexY][11] * offsetY2 * 3 * offsetX2 +
                coefs[indexX][indexY][13] * offsetY3 +
                coefs[indexX][indexY][14] * offsetY3 * 2 * offsetX +
                coefs[indexX][indexY][15] * offsetY3 * 3 * offsetX2;

        s_tag_y = coefs[indexX][indexY][4] +
                coefs[indexX][indexY][5] * offsetX +
                coefs[indexX][indexY][6] * offsetX2 +
                coefs[indexX][indexY][7] * offsetX3 +
                coefs[indexX][indexY][8] * 2 * offsetY +
                coefs[indexX][indexY][9] * 2 * offsetY * offsetX +
                coefs[indexX][indexY][10] * 2 * offsetY * offsetX2 +
                coefs[indexX][indexY][11] * 2 * offsetY * offsetX3 +
                coefs[indexX][indexY][12] * 3 * offsetY2 +
                coefs[indexX][indexY][13] * 3 * offsetY2 * offsetX +
                coefs[indexX][indexY][14] * 3 * offsetY2 * offsetX2 +
                coefs[indexX][indexY][15] * 3 * offsetY2 * offsetX3;
        return true;
    }

    public String toString() {
     return "Spline2D "+breakIntervalX+" "+breaksX.length+" "+breakIntervalY+" "+breaksY.length+" "+s+" "+s_tag_x+" "+s_tag_y;
    }
    private final boolean calcEvenBreaks(double x, double y) {
        int indexX = (int) ((x - breaksX[0]) / breakIntervalX);
        int indexY = (int) ((y - breaksY[0]) / breakIntervalY);
        double extraX,extraY;
        //if ((indexX >= coefs.length) || (indexY >= coefs[0].length)) return false;
        if (indexX >= coefs.length) {
            extraX = x-coefs.length*breakIntervalX;
            indexX = coefs.length-1;
           // meshi.util.Utils.printDebug(this,"X = "+x+" extraX = "+extraX);
        } else extraX = 0;
        if (indexY >= coefs[0].length) {
            extraY = y - coefs[0].length * breakIntervalY;
            indexY = coefs[0].length-1;
          //  meshi.util.Utils.printDebug(this,"Y = "+y+" extraY = "+extraY);
        } else extraY = 0;
        double offsetX = x - breaksX[0] - indexX * breakIntervalX;
        double offsetX2 = offsetX * offsetX;
        double offsetX3 = offsetX2 * offsetX;
        double offsetY = y - breaksY[0] - indexY * breakIntervalY;
        double offsetY2 = offsetY * offsetY;
        double offsetY3 = offsetY2 * offsetY;
//        s = offsetX *
//                (
//                        coefs[indexX][indexY][1] +
//                        coefs[indexX][indexY][2] * offsetX+
//                        coefs[indexX][indexY][3] * offsetX2+
//                        coefs[indexX][indexY][5] * offsetY+
//                        coefs[indexX][indexY][6] * offsetY*offsetX+
//                        coefs[indexX][indexY][7] * offsetY*offsetX2+
//                        coefs[indexX][indexY][9] * offsetY2+
//                        coefs[indexX][indexY][10] * offsetY2*offsetX+
//                        coefs[indexX][indexY][11] * offsetY2*offsetX2+
//                        coefs[indexX][indexY][13] * offsetY3+
//                        coefs[indexX][indexY][14] * offsetY3*offsetX+
//                        coefs[indexX][indexY][15] * offsetY3*offsetX2
//                ) + coefs[indexX][indexY][0] + coefs[indexX][indexY][4] * offsetY + coefs[indexX][indexY][8] * offsetY2 + coefs[indexX][indexY][12] * offsetY3;
//        s = offsetX *
//                (coefs[indexX][indexY][1] +
//                 offsetX * (coefs[indexX][indexY][2] + coefs[indexX][indexY][3] * offsetX) +
//                 offsetY * (coefs[indexX][indexY][5] +
//                            offsetX * (coefs[indexX][indexY][6] +
//                                       coefs[indexX][indexY][7] *  offsetX+
//                                       coefs[indexX][indexY][10] * offsetY+
//                                       coefs[indexX][indexY][11] * offsetY*offsetX+
//                                       coefs[indexX][indexY][14] * offsetY2+
//                                       coefs[indexX][indexY][15] * offsetY2*offsetX) +
//                            offsetY * (coefs[indexX][indexY][9] * offsetY + coefs[indexX][indexY][13] * offsetY))
//                ) + coefs[indexX][indexY][0] +
//                offsetY * (coefs[indexX][indexY][4] +
//                           offsetY * (coefs[indexX][indexY][8] + coefs[indexX][indexY][12] * offsetY))+
//                20*extraX*extraX + 20*extraY*extraY;
        s =     coefs[indexX][indexY][0] +
                coefs[indexX][indexY][1] * offsetX +
                coefs[indexX][indexY][2] * offsetX2 +
                coefs[indexX][indexY][3] * offsetX3 +
                coefs[indexX][indexY][4] * offsetY +
                coefs[indexX][indexY][5] * offsetY * offsetX +
                coefs[indexX][indexY][6] * offsetY * offsetX2 +
                coefs[indexX][indexY][7] * offsetY * offsetX3 +
                coefs[indexX][indexY][8] * offsetY2 +
                coefs[indexX][indexY][9] * offsetY2 * offsetX +
                coefs[indexX][indexY][10] * offsetY2 * offsetX2 +
                coefs[indexX][indexY][11] * offsetY2 * offsetX3 +
                coefs[indexX][indexY][12] * offsetY3 +
                coefs[indexX][indexY][13] * offsetY3 * offsetX +
                coefs[indexX][indexY][14] * offsetY3 * offsetX2 +
                coefs[indexX][indexY][15] * offsetY3 * offsetX3 +
                20*extraX*extraX + 20*extraY*extraY;

        s_tag_x = coefs[indexX][indexY][1] +
                coefs[indexX][indexY][2] * 2 * offsetX +
                coefs[indexX][indexY][3] * 3 * offsetX2 +
                coefs[indexX][indexY][5] * offsetY +
                coefs[indexX][indexY][6] * offsetY * 2 * offsetX +
                coefs[indexX][indexY][7] * offsetY * 3 * offsetX2 +
                coefs[indexX][indexY][9] * offsetY2 +
                coefs[indexX][indexY][10] * offsetY2 * 2 * offsetX +
                coefs[indexX][indexY][11] * offsetY2 * 3 * offsetX2 +
                coefs[indexX][indexY][13] * offsetY3 +
                coefs[indexX][indexY][14] * offsetY3 * 2 * offsetX +
                coefs[indexX][indexY][15] * offsetY3 * 3 * offsetX2 +
                40*extraX;

        s_tag_y = coefs[indexX][indexY][4] +
                coefs[indexX][indexY][5] * offsetX +
                coefs[indexX][indexY][6] * offsetX2 +
                coefs[indexX][indexY][7] * offsetX3 +
                coefs[indexX][indexY][8] * 2 * offsetY +
                coefs[indexX][indexY][9] * 2 * offsetY * offsetX +
                coefs[indexX][indexY][10] * 2 * offsetY * offsetX2 +
                coefs[indexX][indexY][11] * 2 * offsetY * offsetX3 +
                coefs[indexX][indexY][12] * 3 * offsetY2 +
                coefs[indexX][indexY][13] * 3 * offsetY2 * offsetX +
                coefs[indexX][indexY][14] * 3 * offsetY2 * offsetX2 +
                coefs[indexX][indexY][15] * 3 * offsetY2 * offsetX3 +
                40*extraY;
        return true;
    }

    public final String test(double x, double y) {
        double DX = 4e-7;
        double diff;
        calc(x,y);
        double s1 = s;
        double dSdX = s_tag_x;
        double dSdY = s_tag_y;

        x += DX;
        calc(x,y);
        diff = dSdX-(s - s1)/DX;
        if (Math.abs(diff) > testDiff) {
            testDiff = Math.abs(diff);
            return "Testing "+this+" --*-- diff x =  "+diff;
        }

        x -= 2*DX;
        calc(x,y);
        diff = dSdX-(s - s1)/(-DX);
        if (Math.abs(diff) > testDiff) {
            testDiff = Math.abs(diff);
            return "Testing "+this+" --*-- diff -x =  "+diff;
        }

        x += DX;

        y += DX;
        calc(x,y);
        diff = dSdY-(s - s1)/DX;
        if (Math.abs(diff) > testDiff) {
            testDiff = Math.abs(diff);
            return "Testing "+this+" --*-- diff y =  "+diff;
        }

        y -= 2*DX;
        calc(x,y);
        diff = dSdY-(s - s1)/(-DX);
        if (Math.abs(diff) > testDiff) {
            testDiff = Math.abs(diff);
            return "Testing "+this+" --*-- diff -y =  "+diff;
        }

       return null;
    }
    public final void test(MeshiWriter writer, double x, double y) {
        writer.println("testing splin2D with "+x+" , "+y);

        int indexX = (int) ((x - breaksX[0]) / breakIntervalX);
        int indexY = (int) ((y - breaksY[0]) / breakIntervalY);
        double offsetX = x - breaksX[0] - indexX * breakIntervalX;
        double offsetX2 = offsetX * offsetX;
        double offsetX3 = offsetX2 * offsetX;
        double offsetY = y - breaksY[0] - indexY * breakIntervalY;
        double offsetY2 = offsetY * offsetY;
        double offsetY3 = offsetY2 * offsetY;
        writer.println("Offsets = "+indexX+" , "+indexY+" , "+offsetX+" , "+offsetY);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                int ij = 4*i+j;
                writer.println(i+" , "+j+" , "+ij+" , "+coefs[indexX][indexY][ij]);
            }
 //       k   l   4k+l  coef
//        0 , 0 , 0 , -0.3941449169
//        0 , 1 , 1 , 0.7265956463
//        0 , 2 , 2 , -0.1503763176
//        0 , 3 , 3 , 0.0803867053
//        1 , 0 , 4 , -0.9985350248
//        1 , 1 , 5 , -0.468020461
//        1 , 2 , 6 , -0.0944823646
//        1 , 3 , 7 , 0.0999889548
//        2 , 0 , 8 , -0.5617472431
//        2 , 1 , 9 , -1.150938236
//        2 , 2 , 10 , 0.5337753423
//        2 , 3 , 11 , -0.0416629154
//        3 , 0 , 12 , 0.6041198775
//        3 , 1 , 13 , 0.800263096
//        3 , 2 , 14 , -0.2374863278
//        3 , 3 , 15 , -0.0625123305
        s =     coefs[indexX][indexY][0] +                       // x^0 y^0
                coefs[indexX][indexY][1] * offsetX +             // x^1 y^0
                coefs[indexX][indexY][2] * offsetX2 +            // x^2 y^0
                coefs[indexX][indexY][3] * offsetX3 +            // x^3 y^0
                coefs[indexX][indexY][4] * offsetY +             // x^0 y^1
                coefs[indexX][indexY][5] * offsetY * offsetX +   // x^0 y^1
                coefs[indexX][indexY][6] * offsetY * offsetX2 +
                coefs[indexX][indexY][7] * offsetY * offsetX3 +
                coefs[indexX][indexY][8] * offsetY2 +
                coefs[indexX][indexY][9] * offsetY2 * offsetX +
                coefs[indexX][indexY][10] * offsetY2 * offsetX2 +
                coefs[indexX][indexY][11] * offsetY2 * offsetX3 +
                coefs[indexX][indexY][12] * offsetY3 +
                coefs[indexX][indexY][13] * offsetY3 * offsetX +
                coefs[indexX][indexY][14] * offsetY3 * offsetX2 +
                coefs[indexX][indexY][15] * offsetY3 * offsetX3;

        s_tag_x = coefs[indexX][indexY][1] +
                coefs[indexX][indexY][2] * 2 * offsetX +
                coefs[indexX][indexY][3] * 3 * offsetX2 +
                coefs[indexX][indexY][5] * offsetY +
                coefs[indexX][indexY][6] * offsetY * 2 * offsetX +
                coefs[indexX][indexY][7] * offsetY * 3 * offsetX2 +
                coefs[indexX][indexY][9] * offsetY2 +
                coefs[indexX][indexY][10] * offsetY2 * 2 * offsetX +
                coefs[indexX][indexY][11] * offsetY2 * 3 * offsetX2 +
                coefs[indexX][indexY][13] * offsetY3 +
                coefs[indexX][indexY][14] * offsetY3 * 2 * offsetX +
                coefs[indexX][indexY][15] * offsetY3 * 3 * offsetX2;

        s_tag_y = coefs[indexX][indexY][4] +
                coefs[indexX][indexY][5] * offsetX +
                coefs[indexX][indexY][6] * offsetX2 +
                coefs[indexX][indexY][7] * offsetX3 +
                coefs[indexX][indexY][8] * 2 * offsetY +
                coefs[indexX][indexY][9] * 2 * offsetY * offsetX +
                coefs[indexX][indexY][10] * 2 * offsetY * offsetX2 +
                coefs[indexX][indexY][11] * 2 * offsetY * offsetX3 +
                coefs[indexX][indexY][12] * 3 * offsetY2 +
                coefs[indexX][indexY][13] * 3 * offsetY2 * offsetX +
                coefs[indexX][indexY][14] * 3 * offsetY2 * offsetX2 +
                coefs[indexX][indexY][15] * 3 * offsetY2 * offsetX3;
    }

    private boolean testForEvenBreaks(double[] vec) {
        double interval = vec[1] - vec[0];
        for (int c = 1; c < vec.length; c++) {//            meshi.util.Utils.printDebug(this,"interval = "+interval+", "+((vec[c] - vec[0]) - interval)+", "+(interval + (vec[c] -vec[0])));
//            if (Math.abs((vec[c] - vec[0]) - interval) / (interval + (vec[c] -vec[0])) > 1e-13) {
            if (Math.abs((vec[c] - vec[c-1]) - interval) / (interval + (vec[c] -vec[c-1])) > 1e-13) {
               return false;
            }
        }
		return true;
	}
}	
