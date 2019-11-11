package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore;

public class BoreCalculator {
private  final double K = 1;
private  double MINx, MAXx, slope, leftLineBound, rightLineBound;

public static double f, df;

 public BoreCalculator(double MINx, double MAXx, double slope) {
    this.MINx = MINx;
    this.MAXx = MAXx;
    this.slope = slope;    
    leftLineBound = MINx+slope;
    rightLineBound = MAXx-slope;
 }

    /**
    * energy and dirivarives calculation.
    **/
  protected  void calculateBore(double x) {
      double d, d2;
      double a, b;

      if ((x < MINx) || (x >  MAXx)) {
          f = 0;
          df = 0;
          return;
      };
      if (x < leftLineBound)  {
          a = -2./(slope*slope*slope);
          b = 3./(slope*slope);
          d = x - MINx;
          d2 = d*d;
          f = K*(a*d2*d + b*d2);
          df = K*(3*a*d2 + 2*b*d);
      }
      else
      if (x > rightLineBound) {
          a = 2./(slope*slope*slope);
          b = -3./(slope*slope);
          d = x-rightLineBound;
          d2 = d*d;
          f = K*(a*d2*d + b*d2 + 1);
          df = K*(3*a*d2 + 2*b*d);
      }
      else {
          f = K;
          df = 0;
      }
  }


 protected double df(){return df;}
 protected double f(){return f;}
   
 public String toString() {
       return "BoreElement: energy = " + f + "derivative = " + df;
   }
 
}
