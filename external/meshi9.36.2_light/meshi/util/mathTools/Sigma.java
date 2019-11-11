/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.mathTools;

import meshi.util.*;

/**
 * This class host the method "sigma" which gives a derivable sigmoid function
 * of the following form:
 * <p/>
 * x<0   - not defined
 * 0<x<p1 - a linear descent from 1.0 to the value of valAtp1
 * p1<x<p2 - a cubic spline descent from (p1,valAtp1) to (p2,valAtp2)
 * p2<x<end - a quadratic descent from (p2,valAtp2) to (end,0.0).
 * end<x - 0.0
 * <p/>
 * After the method sigma(...) is run, the public fields s and s_tag are updated:
 * s - the sigmoid value at x.
 * s_tag - the sigmoid first derivative at x.
 * <p/>
 * Note: The sigmoid function is derivable for all x=<0.
 */

public class Sigma {
    static final double EPSILON = 0.2;
    private double s; // The value of sigma
    private double s_tag; // The derivative of sigma
    public  final double s() {return s;}
    public  final double s_tag() {return s_tag;}
    public double end, p1, p2, valAtp1, valAtp2;

    public Sigma() {
        s = 0;
        s_tag = 0;
    }

    public String toString(){
        return "Sigma: "+end+" "+p1+" "+p2+" "+valAtp1+" "+valAtp2;
    }

    public Sigma(double end, double p1, double p2,
                 double valAtp1, double valAtp2) {
        this.end = end;
        this.p1 = p1;
        this.p2 = p2;
        this.valAtp1 = valAtp1;
        this.valAtp2 = valAtp2;
        calcABCD();
    }

        double a1, a2, v1, v2, a, b, c, d; //Auxilary variables

    public void calcABCD() {
        v1 = valAtp1;
        v2 = valAtp2;
        a1 = (valAtp1 - 1) / p1;
        a2 = 2 * valAtp2 / (p2 - end);
        a = (p1 * a1 - 2 * v1 - p2 * a1 - p2 * a2 + p1 * a2 + 2 * v2) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
        b = -(2 * p1 * p1 * a2 + p1 * p1 * a1 + 3 * p1 * v2 - p1 * p2 * a2 - 3 * p1 * v1 + p1 * p2 * a1 - 2 * p2 * p2 * a1 + 3 * p2 * v2 - p2 * p2 * a2 - 3 * p2 * v1) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
        c = (2 * p1 * p1 * p2 * a1 + p1 * p1 * p2 * a2 + p1 * p1 * p1 * a2 - p2 * p2 * p1 * a1 + 6 * p1 * p2 * v2 - 2 * p2 * p2 * p1 * a2 - 6 * p1 * p2 * v1 - p2 * p2 * p2 * a1) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
        d = -(p1 * p1 * p1 * p2 * a2 - p1 * p1 * p1 * v2 + p1 * p1 * p2 * p2 * a1 + 3 * p1 * p1 * p2 * v2 - p1 * p1 * p2 * p2 * a2 - p1 * a1 * p2 * p2 * p2 + v1 * p2 * p2 * p2 - 3 * v1 * p1 * p2 * p2) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
    }


    //sigma.sigma(1.0 + cosAng, 2, 1.0 + angSigmoidBegins, 1.0 + angSigmoidEnds, 1.0, 0.0);
    //                 x        end  p1                       p2                 val@p1 val@p2
    // in 180          0        2    0.83                   1                    1     0
    public final void sigma(double x,double end, double p1, double p2,
                            double valAtp1, double valAtp2) {
        this.end = end;
        this.p1 = p1;
        this.p2 = p2;
        this.valAtp1 = valAtp1;
        this.valAtp2 = valAtp2;
        calcABCD();
        sigma(x);
    }
    public final void sigma(double x) {


        if (x >= end) {
            s = 0;
            s_tag = 0;
        }
        else if (x >= p2) {
                s = (x - end) * (x - end) * valAtp2 / ((p2 - end) * (p2 - end));
                s_tag = 2 * (x - end) * valAtp2 / ((p2 - end) * (p2 - end));
        }
        else if (x >= p1) {
            double x2 = x*x;
            double ax2 = a*x2;
            double bx = b*x;
            //calcABCD(end, p1, p2, valAtp1, valAtp2);
            s = ax2 * x + bx * x + c * x + d;
            s_tag = 3 * ax2 + 2 * bx + c;

        }
//                v1 = valAtp1;
//                v2 = valAtp2;
//                a1 = (valAtp1 - 1) / p1;
//                a2 = 2 * valAtp2 / (p2 - end);
//            double a1PlusA2 = a1 + a2;
//            double p1PlusP2 = p1 + p2;
//            double p1MinusP2 = p1 - p2;
//            double v1MinusV2 = v1 - v2;
//            double p1p2      = p1*p2;
//            double p2p2      = p2*p2;
//            double p1p1      = p1*p1;
//            double p2p2p2    = p2p2*p2;
//            double p1p1p1    = p1p1*p1;
//            double Q         = -3 *p1p2*p1MinusP2 - p2p2p2 + p1p1p1;
//
//              //a = (p1 * a1 - 2 * v1 - p2 * a1 - p2 * a2 + p1 * a2 + 2 * v2) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
////              a = (p1 * (a1 + a2) - p2 *(a1 + a2) - 2 *( v1   - v2)) / (3 * p1 * p2 * (p2 - p1) - p2 * p2 * p2 + p1 * p1 * p1  );
////                a = ((a1 + a2) * (p1 - p2) - 2 *( v1 - v2)) / (3 * p1 * p2 * (p2 - p1) - p2 * p2 * p2 + p1 * p1 * p1  );
//            a = (a1PlusA2*p1MinusP2 - 2*v1MinusV2) / Q;
////            if (a != aNew) {
////            }
//              //b = -(2 * p1 * p1 * a2 + p1 * p1 * a1 + 3 * p1 * v2 - p1 * p2 * a2 - 3 * p1 * v1 + p1 * p2 * a1 - 2 * p2 * p2 * a1 + 3 * p2 * v2 - p2 * p2 * a2 - 3 * p2 * v1) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
////                double bNew = -(p1 * p1 * (2 * a2 + a1)  + p1 * p2 * (a1 - a2) - p2 * p2 * (2 * a1 + a2) + 3 * p1 * (v2  - v1)  + 3 * p2 * (v2  - v1)) /
////                               ( 3 * p1 * p2 * (p2 - p1) - p2 * p2 * p2  + p1 * p1 * p1);
//             b = -(p1p1 * (2 * a2 + a1)  + p1p2 * (a1 - a2) - p2p2 * (2 * a1 + a2) - 3 * p1PlusP2 * v1MinusV2  ) / Q;
////            if (b != bNew) {
////            }
//                c = (2 * p1 * p1 * p2 * a1 + p1 * p1 * p2 * a2 + p1 * p1 * p1 * a2 - p2 * p2 * p1 * a1 + 6 * p1 * p2 * v2 - 2 * p2 * p2 * p1 * a2 - 6 * p1 * p2 * v1 - p2 * p2 * p2 * a1) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
//            double cNew  = (p1p1 * p2 * (2 * a1 + a2) + p1p1p1 * a2 - p2p2 * p1 * (2 * a2 + a1) - 6 * p1p2 * v1MinusV2 - p2p2p2 * a1) / Q;
////            if (c != cNew) {
////            }
//
//            d = -(p1 * p1 * p1 * p2 * a2 - p1 * p1 * p1 * v2 + p1 * p1 * p2 * p2 * a1 + 3 * p1 * p1 * p2 * v2 - p1 * p1 * p2 * p2 * a2 - p1 * a1 * p2 * p2 * p2 + v1 * p2 * p2 * p2 - 3 * v1 * p1 * p2 * p2) / (-p2 * p2 * p2 - 3 * p2 * p1 * p1 + p1 * p1 * p1 + 3 * p1 * p2 * p2);
//                s = a * x * x * x + b * x * x + c * x + d;
//                s_tag = 3 * a * x * x + 2 * b * x + c;
//        }
        else if (x >= 0.0) {
            s = 1 - (1 - valAtp1) / p1 * x;
            s_tag = -(1 - valAtp1) / p1;
        }
        else throw new RuntimeException("A negative x-value is not possible. Yet x = " + x);
        //chen 3.10.12
        //flattening the lower part of the sigma function to reduce the step
        //height during Hb-formation.
        double epsilonPlus1 = EPSILON + 1;
        double epsilonPlusS = EPSILON + s;
        s_tag = s_tag*s*epsilonPlus1*(2 - s/epsilonPlusS)/epsilonPlusS;
        s     = epsilonPlus1*s*s/epsilonPlusS;
    }

}