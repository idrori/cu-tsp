/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.Atom;

public class ViewAt {
    public ViewAt() {
    }

    public static double[][] transformToOrigin(double M[][], double invM[][], Atom P1, Atom P2, Atom P3) {
        double[] p1 = new double[4];
        double[] p2 = new double[4];
        double[] p3 = new double[4];

        p1[0] = P1.x();
        p1[1] = P1.y();
        p1[2] = P1.z();

        p2[0] = P2.x();
        p2[1] = P2.y();
        p2[2] = P2.z();

        p3[0] = P3.x();
        p3[1] = P3.y();
        p3[2] = P3.z();

        p1[3] = p2[3] = p3[3] = 1.0;
        return transformToOrigin(M, invM, p1, p2, p3);
    }

    public static double[][] transformToOrigin(double M[][], double invM[][], double P1[], double P2[], double P3[]) {
        // double invM[][] = new double[4][4];
        double d12;
        double P120, P121, P122, P130, P131, P132;
        P120 = P2[0] - P1[0];
        P121 = P2[1] - P1[1];
        P122 = P2[2] - P1[2];
        P130 = P3[0] - P1[0];
        P131 = P3[1] - P1[1];
        P132 = P3[2] - P1[2];
        d12 = Math.sqrt(P120 * P120 + P121 * P121 + P122 * P122);
        invM[2][0] = M[0][2] = P120 / d12;
        invM[2][1] = M[1][2] = P121 / d12;
        invM[2][2] = M[2][2] = P122 / d12;

        double[] tmp = nmcrosprod(P130, P131, P132, P120, P121, P122);
        M[0][0] = tmp[0];
        M[1][0] = tmp[1];
        M[2][0] = tmp[2];

        invM[0][0] = M[0][0];
        invM[0][1] = M[1][0];
        invM[0][2] = M[2][0];

        M[0][3] = invM[0][3] = 0;
        M[1][3] = invM[1][3] = 0;
        M[2][3] = invM[2][3] = 0;

        tmp = nmcrosprod(M[0][2], M[1][2], M[2][2], M[0][0], M[1][0], M[2][0]);
        M[0][1] = tmp[0];
        M[1][1] = tmp[1];
        M[2][1] = tmp[2];

        invM[1][0] = M[0][1];
        invM[1][1] = M[1][1];
        invM[1][2] = M[2][1];

        invM[3][0] = P1[0];
        invM[3][1] = P1[1];
        invM[3][2] = P1[2];
        invM[3][3] = 1.0;
//Matrix M = new Matrix(invM);
        M[3][0] = -P1[0] * M[0][0] - P1[1] * M[1][0] - P1[2] * M[2][0];
        M[3][1] = -P1[0] * M[0][1] - P1[1] * M[1][1] - P1[2] * M[2][1];
        M[3][2] = -P1[0] * M[0][2] - P1[1] * M[1][2] - P1[2] * M[2][2];
        M[3][3] = 1.0;
/*
  System.out.println("M:");
  for (int i=0;i<4;i++){
      for (int j=0;j<4;j++)
          System.out.print(M[i][j]+" ");
      System.out.println();
  }
  System.out.println("invM:");
  for (int i=0;i<4;i++){
      for (int j=0;j<4;j++)
          System.out.print(invM[i][j]+" ");
      System.out.println();
  }
  System.out.println();
         double	M[4][4], invM[4][4];
                  double P1[4], P2[4], P3[4]; lookfrom,lookat,lookup

 */
        return invM;
    }

    protected static double[] nmcrosprod(double x1, double y1, double z1, double x2, double y2, double z2) {
/*     double	x1, y1, z1, x2, y2, z2; r1 cross r2
     double	*x3, *y3, *z3;		 Normalized crossproduct
*/
        double dis;        /* length of crossproduct vector r1 x r2 */
        double[] xyz = new double[3];
        double x, y, z;

        x = y1 * z2 - y2 * z1;
        y = z1 * x2 - z2 * x1;
        z = x1 * y2 - x2 * y1;

        dis = Math.sqrt((x * x + y * y + z * z));

        xyz[0] = (x / dis);
        xyz[1] = (y / dis);
        xyz[2] = (z / dis);

        return xyz;
    }


    /*protected  double  dis (double x,double y,double z){
          return (Math.sqrt( (x*x + y*y + z*z) ) )  ;
      }
    */
/*
  *x3 = x / dis;
  *y3 = y / dis;
  *z3 = z / dis;  */


}
