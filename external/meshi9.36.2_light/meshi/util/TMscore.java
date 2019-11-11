package meshi.util;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.sequences.*;

import javax.swing.*;
import javax.swing.text.html.HTMLDocument;
import java.util.Iterator;

/**
 * Created by chen on 05/03/2017.
 */
public class TMscore {
    static int nmax = 3000;
    static int n_cut;// ![1,n_ali], align residues for the score
    double[] xt = new double[nmax];
    double[] yt = new double[nmax];
    double[] zt = new double[nmax];
    double[] w = new double[nmax];
    double[][] u = new double[4][4];
    double[] t = new double[4];
    static int[] i_ali = new int[nmax];
    static int[] k_ali = new int[nmax];
    static int[] k_ali0 = new int[nmax];
    int[] L_ini = new int[100];
    int[] iq = new int[nmax];
    int n_it = 20;                   //!maximum number of iterations
    int d_output = 5;                //!for output alignment
    int n_init_max = 6;               //!maximum number of L_init
    int n_init = 0;
    int L_ini_min = 4;
    int L_init;
    int iL_max;
    static int n_ali;
    int LL;
    int ka, ka0;
    double[][] r_1 = new double[4][nmax + 1];
    double[][] r_2 = new double[4][nmax + 1];
    static double[] xa = new double[nmax];
    static double[] ya = new double[nmax];
    static double[] za = new double[nmax];
    static double[] xb = new double[nmax];
    static double[] yb = new double[nmax];
    static double[] zb = new double[nmax];
    static int[] nresA, nresB;
    static int[] iA = new int[nmax];
    static int[] iB = new int[nmax];;
    static int nseqA;
    static int nseqB;
    static double d0, d0_search;

    public TMscore(Protein model, Protein nativeStructure) throws AlignmentException {

        int L1 = model.chain().numberOfNonDummyResidues();
        nseqA = L1;
        xa = new double[L1 + 1];
        ya = new double[L1 + 1];
        za = new double[L1 + 1];
        nresA = new int[L1 + 1];
        Iterator resIter = model.chain().iterator();
        Residue residue = null;
        int[] residueIndexA = new int[model.chain().size()];
        for (int i = 1; i <= L1; i++) {
            while (resIter.hasNext() && (residue = (Residue) resIter.next()).dummy()) ;
            if (residue.dummy())
                throw new RuntimeException("This is weird");
            Atom ca = residue.ca();
            xa[i] = ca.x();
            ya[i] = ca.y();
            za[i] = ca.z();
            nresA[i] = residue.number();
            residueIndexA[residue.number()] = i;
        }


        int L2 = nativeStructure.chain().numberOfNonDummyResidues();
        nseqB = L1;
        xb = new double[L2 + 1];
        yb = new double[L2 + 1];
        zb = new double[L2 + 1];
        int[] nresB = new int[L2 + 1];
        resIter = nativeStructure.chain().iterator();
        residue = null;
        int[] residueIndexB = new int[nativeStructure.chain().size()];

        for (int i = 1; i <= L2; i++) {
            while (resIter.hasNext() && (residue = (Residue) resIter.next()).dummy()) ;
            Atom ca = residue.ca();
            xb[i] = ca.x();
            yb[i] = ca.y();
            zb[i] = ca.z();
            nresB[i] = residue.number();
            residueIndexB[residue.number()] = i;
        }

        ResidueAlignment alignment = new ResidueAlignment(model.chain(), model.name(), nativeStructure.chain(), nativeStructure.name(), ResidueAlignmentMethod.IDENTITY);
        Utils.println("Calculating TM-score with alignment \n"+alignment.getSequenceAlignment().toString()+"\n"+alignment.toString());
        n_ali = alignment.size();
        //Utils.printDebug(this,"zzzzzz "+n_ali);
        //int[] iA = new int[n_ali];
        //int[] iB = new int[n_ali];


        for (int i = 1; i < n_ali; i++) {
            ResidueAlignmentColumn column = alignment.get(i);
            Residue modelResidue = ((ResidueAlignmentCell) column.cell0()).residue();
            Residue nativeResidue = ((ResidueAlignmentCell) column.cell1()).residue();
            iA[i] = residueIndexA[modelResidue.number()];
            iB[i] = residueIndexB[nativeResidue.number()];
        }
        /************/////
        /*     parameters:
        *****************/
        /***   d0-------------> */
        d0 = 1.24 * Math.pow(nseqB - 15.0, 1.0 / 3.0) - 1.8;  //d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
        if (d0 < 0.5) d0 = 0.5;                            //if(d0.lt.0.5)d0=0.5
        //***   d0_search ----->
        d0_search = d0;
        if (d0_search > 8) d0_search = 8; //if(d0_search.gt.8)d0_search=8
        if (d0_search < 4.5) d0_search = 4.5; //if(d0_search.lt.4.5)d0_search=4.5
        //***   iterative parameters ----->
        if (n_ali < 4) L_ini_min = n_ali; //if(n_ali.lt.4)L_ini_min=n_ali
        for (int i = 1; i <= n_init_max - 1; i++) { // do i=1,n_init_max-1
            n_init = n_init + 1;
            L_ini[n_init] = (int) ((1.0 * n_ali) / Math.pow(2.0, (1.0 * n_init - 1.0))); // L_ini(n_init) = n_ali / 2 **(n_init - 1)
            if (L_ini[n_init] <= L_ini_min) { //if (L_ini(n_init).le.L_ini_min) then
                L_ini[n_init] = L_ini_min;
                break; //goto 402
            }//endif
        }// enddo
        n_init = n_init + 1;
        L_ini[n_init] = L_ini_min;

        //402  continue

    }

    public double findMaxScore() {
        /******************************************************************
         *     find the maximum score starting from local structures superposition
         ********************************************************************/
        double score_max = -1;              //TM-score
        for (int i_init = 1; i_init <= n_init; i_init++) { //do 333 i_init=1,n_init
            L_init = L_ini[i_init];
            iL_max = n_ali - L_init + 1;
            for (int iL = 1; iL <= iL_max; iL++) { //do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
                LL = 0;
                ka = 0;
               // Utils.printDebug(this,"xxxxxxxxxx "+i_init+" "+iL_max+" "+L_init+" "+iL+" "+(iL + L_init - 1));
                for (int i = 1; i <= L_init; i++) {//do i=1,L_init
                    int k = iL + i - 1;    //          ![1,n_ali] common aligned
                    r_1[1][i] = xa[iA[k]];
                    r_1[2][i] = ya[iA[k]];
                    r_1[3][i] = za[iA[k]];
                    r_2[1][i] = xb[iB[k]];
                    r_2[2][i] = yb[iB[k]];
                    r_2[3][i] = zb[iB[k]];
                    LL = LL + 1;
                    ka = ka + 1;
                    k_ali[ka] = k;
                } //enddo
                U3b u3b = new U3b(r_1, r_2, LL); //!u rotate r_1 to r_2// call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
                u3b.calc();
                double armsd;
                double Rcomm;
                //if (i_init == 1) { //!global superposition
                    armsd = Math.sqrt(u3b.rms / LL);
                    Rcomm = armsd;
               // }
                for (int j = 1; j <= nseqA; j++) {//do j = 1, nseqA
                    xt[j] = u3b.t[1] + u3b.u[1][1] * xa[j] + u3b.u[1][2] * ya[j] + u3b.u[1][3] * za[j];
                    yt[j] = u3b.t[2] + u3b.u[2][1] * xa[j] + u3b.u[2][2] * ya[j] + u3b.u[2][3] * za[j];
                    zt[j] = u3b.t[3] + u3b.u[3][1] * xa[j] + u3b.u[3][2] * ya[j] + u3b.u[3][3] * za[j];
                } //enddo
                for (int j = 1; j <= nseqA; j++) {//do j = 1, nseqA
                    xa[j] = xt[j];
                    ya[j] = yt[j];
                    za[j] = zt[j];
                }
                    double d = d0_search - 1;
                ScoreFun scoreFun = new ScoreFun(armsd < 2);//call score_fun !init, get scores, n_cut + i_ali[i] for iteration
                double score = scoreFun.score;
                //Utils.printDebug(this,"yyyyyyyyyyyyyyyyyyyyy "+armsd+" "+LL+" "+score+" "+score_max);
                if (score_max < score) {//then
                    score_max = score;
                    ka0 = ka;
                    for (int i = 1; i <= ka0; i++) {//do i = 1, ka0
                        k_ali0[i] = k_ali[i];
                    }//enddo
                }//endif
                /***iteration for extending---------------------------------- > */
                d = d0_search + 1;
                for (int it = 1; it <= n_it; it++) {//do 301 it = 1, n_it
                    LL = 0;
                    ka = 0;
                    for (int i = 1; i <= n_cut; i++) {//do i = 1, n_cut
                        int m = i_ali[i];
                        ;// ![1][n_ali];
                        r_1[1][i] = xa[iA[m]];
                        r_1[2][i] = ya[iA[m]];
                        r_1[3][i] = za[iA[m]];
                        r_2[1][i] = xb[iB[m]];
                        r_2[2][i] = yb[iB[m]];
                        r_2[3][i] = zb[iB[m]];
                        ka = ka + 1;
                        k_ali[ka] = m;
                        LL = LL + 1;
                    }//enddo
                    u3b = new U3b(r_1, r_2, LL);//call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
                    u3b.calc();
                    armsd = Math.sqrt(u3b.rms / LL);
                    for (int j = 1; j <= nseqA; j++) {//do j = 1, nseqA
                        xt[j] = u3b.t[1] + u3b.u[1][1] * xa[j] + u3b.u[1][2] * ya[j] + u3b.u[1][3] * za[j];
                        yt[j] = u3b.t[2] + u3b.u[2][1] * xa[j] + u3b.u[2][2] * ya[j] + u3b.u[2][3] * za[j];
                        zt[j] = u3b.t[3] + u3b.u[3][1] * xa[j] + u3b.u[3][2] * ya[j] + u3b.u[3][3] * za[j];
                    }//enddo
                    for (int j = 1; j <= nseqA; j++) {//do j = 1, nseqA
                        xa[j] = xt[j];
                        ya[j] = yt[j];
                        za[j] = zt[j];
                    }
                    scoreFun = new ScoreFun(armsd < 2);
                    score = scoreFun.score; //call score_fun !get scores, n_cut + i_ali[i] for iteration
                    //Utils.printDebug(this,"yyyyyyyyyyyyyyyyyyyyy2 "+armsd+" "+LL+" "+score+" "+score_max);
                    if (score_max < score) {//then
                        score_max = score;
                        ka0 = ka;
                        for (int i = 1; i <= 3; i++) {//do i = 1, ka
                            k_ali0[i] = k_ali[i];
                        }//enddo
                    }//endif
                    if (it == n_it) break; //goto 302
                    int neq;
                    if (n_cut == ka) {//then
                        neq = 0;
                        for (int i = 1; i <= 3; i++) {//do i = 1, n_cut
                            if (i_ali[i] == k_ali[i]) neq = neq + 1;
                        }//enddo
                        if (n_cut == neq) break; //goto 302
                    }//endif
                }//301 continue ! for iteration
                //302 continue
                //300 continue ! for shift
                //333 continue ! for initial length, L_ali/M
                //return score_max;
            }
        }
//        if (1 == 1)
//            throw new RuntimeException("Debug score_max = "+score_max);
        return score_max;
    }


//---------------------------------------------

    /*cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
    c  w    - w(m) is weight for atom pair   c m                (given)
    c  x    - x(i,m) are coordinates of atom c m in set x       (given)
    c  y    - y(i,m) are coordinates of atom c m in set y       (given)
    c  n    - n is number of atom pairs                         (given)
    c  mode  - 0:calculate rms only                             (given)
    c          1:calculate rms,u,t                              (takes longer)
    c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
    c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
    c  t    - t(i)   is translation vector for best superposition  (result)
    c  ier  - 0: a unique optimal superposition has been determined(result)
    c       -1: superposition is not unique but optimal
    c       -2: no result obtained because of negative weights w
    c           or all weights equal to zero.
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  */
    private static class U3b {
        //double precision w(*), x(3,*), y(3,*)
        //integer n, mode

        //double precision rms, u(3,3), t(3)
        //integer ier
        double[] w = new double[nmax];
        double[][] x = new double[4][nmax];
        double[][] y = new double[4][nmax];
        int n;
        int mode = 1;
        double rms;
        double[][] u = new double[4][4];
        double[] t = new double[4];
        int ier;
        private int i, j, k, l, m1, m;
        double[][] r = new double[4][4];
        double[] xc = new double[4];
        double[] yc = new double[4];
        double wc;
        double[][] a = new double[4][4];
        double[][] b = new double[4][4];
        double[] e = new double[4];
        double[] rr = new double[7];
        double[] ss = new double[7];
        double e0, d, spur, det, cof, h, g;
        double cth, sth, sqrth, p, sigma;


        double sqrt3 = 1.73205080756888d + 00;
        double tol = 1.0d - 2;
        double zero = 0.0d + 00;
        int[] ip = {0, 1, 2, 4, 2, 3, 5, 4, 5, 6};
        int[] ip2312 = {0, 2, 3, 1, 2};

        public U3b(double[][] x, double[][] y, int n) {
            //subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
            //call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
            this.x = x;
            this.y = y;
            this.n = n;
            for (int i = 0; i < w.length; i++) w[i] = 1;
        }

        public void calc() {
            wc = zero;
            rms = zero;
            e0 = zero;

            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                xc[i] = zero;
                yc[i] = zero;
                t[i] = zero;
                for (int j = 1; j <= 3; j++) {//do j = 1, 3
                    r[i][j] = zero;
                    u[i][j] = zero;
                    a[i][j] = zero;
                    if (i == j) {//if (i.eq.j) then
                        u[i][j] = 1.0;
                        a[i][j] = 1.0;
                    } //end if
                } //end do
            } //end do

            ier = -1;
            if (n < 1) return;//if (n.lt. 1 )return
            ier = -2;
            for (int m = 1; m <= n; m++) {//do m = 1, n
                if (w[m] < 0) return; //if (w(m).lt. 0.0 )return
                wc = wc + w[m];
                for (i = 1; i <= 3; i++) { //do i = 1, 3
                    xc[i] = xc[i] + w[m] * x[i][m];
                    yc[i] = yc[i] + w[m] * y[i][m];
                } //end do
            }//end do
            if (wc <= zero) return;
            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                xc[i] = xc[i] / wc;
                yc[i] = yc[i] / wc;
            } //end do

            for (int m = 1; m <= n; m++) {  //do m = 1, n
                for (int i = 1; i <= 3; i++) {  //do i = 1, 3
                    e0 = e0 + w[m] * (Math.pow((x[i][m] - xc[i]), 2.0) + Math.pow((y[i][m] - yc[i]), 2));
                    d = w[m] * (y[i][m] - yc[i]);
                    for (int j = 1; j <= 3; j++) {  //do j = 1, 3
                        r[i][j] = r[i][j] + d * (x[j][m] - xc[j]);
                    } //end do
                } //end do
            }  //end do

            det = r[1][1] * ((r[2][2] * r[3][3]) - (r[2][3] * r[3][2])) +
                    -r[1][2] * ((r[2][1] * r[3][3]) - (r[2][3] * r[3][1])) +
                    +r[1][3] * ((r[2][1] * r[3][2]) - (r[2][2] * r[3][1]));

            sigma = det;

            m = 0;
            for (int j = 1; j <= 3; j++) {  //do i = 1, 3
                for (int i = 1; i <= j; i++) {  //do i = 1, j
                    m = m + 1;
                    rr[m] = r[1][i] * r[1][j] + r[2][i] * r[2][j] + r[3][i] * r[3][j];
                } //end do
            } //end do

            spur = (rr[1] + rr[3] + rr[6]) / 3.0;
            cof = (((((rr[3] * rr[6] - rr[5] * rr[5]) + rr[1] * rr[6]) +
                    -rr[4] * rr[4]) + rr[1] * rr[3]) - rr[2] * rr[2]) / 3.0;
            det = det * det;
            for (int i = 1; i <= 3; i++) {  //do i = 1, 3
                e[i] = spur;
            }  //end do
            if (spur > zero) {   // if (spur <= zero) goto 40
                d = spur * spur;
                h = d - cof;
                g = (spur * cof - det) / 2.0 - spur * h;
                /*if (h <= zero) { //then
                    if (mode == 0) {//then
                        goto 50
                        else
                    goto 30
                    }//end if
                }end if */
                if (h > zero) {  //30 & 50
                    sqrth = Math.sqrt(h);
                    d = h * h * h - g * g;
                    if (d < zero) d = zero;
                    d = Math.atan2(Math.sqrt(d), -g) / 3.0;
                    cth = sqrth * Math.cos(d);
                    sth = sqrth * sqrt3 * Math.sin(d);
                    e[1] = (spur + cth) + cth;
                    e[2] = (spur - cth) + sth;
                    e[3] = (spur - cth) - sth;

                    /**if (mode ==  0 )then
                     goto 50
                     end     if */
                    if (mode != 0) {
                        int[] lArray = {1, 3, 2};//do l = 1, 3, 2
                        for (int l : lArray) {
                            d = e[l];
                            ss[1] = (d - rr[3]) * (d - rr[6]) - rr[5] * rr[5];
                            ss[2] = (d - rr[6]) * rr[2] + rr[4] * rr[5];
                            ss[3] = (d - rr[1]) * (d - rr[6]) - rr[4] * rr[4];
                            ss[4] = (d - rr[3]) * rr[4] + rr[2] * rr[5];
                            ss[5] = (d - rr[1]) * rr[5] + rr[2] * rr[4];
                            ss[6] = (d - rr[1]) * (d - rr[3]) - rr[2] * rr[2];

                        /*if( dabs(ss(1)) .ge. dabs(ss(3)) ) then
                                j=1
                        if( dabs(ss(1)) .lt. dabs(ss(6)) ) j = 3
                        else if( dabs(ss(3)) .ge. dabs(ss(6)) ) then
                                j = 2
                        else
                            j = 3
                        end if    */

                            if (Math.abs(ss[1]) >= Math.abs(ss[3])) {//then
                                j = 1;
                                if (Math.abs(ss[1]) < Math.abs(ss[6])) j = 3;
                            } else {
                                if (Math.abs(ss[3]) >= Math.abs(ss[6])) //then
                                    j = 2;
                                else
                                    j = 3;
                            }//end if

                            d = zero;
                            j = 3 * (j - 1);

                            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                                k = ip[i + j];
                                a[i][l] = ss[k];
                                d = d + ss[k] * ss[k];
                            }//end do
                            if (d > zero) d = 1.0 / Math.sqrt(d);
                            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                                a[i][l] = a[i][l] * d;
                            }//end do
                        }//end do

                        d = a[1][1] * a[1][3] + a[2][1] * a[2][3] + a[3][1] * a[3][3];
                        if ((e[1] - e[2]) > (e[2] - e[3])) {//then
                            m1 = 3;
                            m = 1;
                        } else {
                            m1 = 1;
                            m = 3;
                        }//endif

                        p = zero;
                        for (int i = 1; i <= 3; i++) {//do i = 1, 3
                            a[i][m1] = a[i][m1] - d * a[i][m];
                            p = p + Math.pow(a[i][m1], 2);
                        }//end do
                        if (p <= tol) {//then
                            p = 1.0;
                            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                                if (p >= Math.abs(a[i][m])) {//if (p < Math.abs(a[i][m])) cycle
                                    p = Math.abs(a[i][m]);
                                    j = i;
                                }
                            }//end do
                            k = ip2312[j];
                            l = ip2312[j + 1];
                            p = Math.sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);//p = Math.sqrt(a[k, m] * * 2 + a[l][m] **2 )
                            if (p > tol) {//if (p <=  tol) goto 40
                                a[j][m1] = zero;
                                a[k][m1] = -a[l][m] / p;
                                a[l][m1] = a[k][m] / p;
                                a[1][2] = a[2][3] * a[3][1] - a[2][1] * a[3][3];
                                a[2][2] = a[3][3] * a[1][1] - a[3][1] * a[1][3];
                                a[3][2] = a[1][3] * a[2][1] - a[1][1] * a[2][3];
                            }
                        } else {
                            p = 1.0 / Math.sqrt(p);
                            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                                a[i][m1] = a[i][m1] * p;
                            }//end do
                            a[1][2] = a[2][3] * a[3][1] - a[2][1] * a[3][3];
                            a[2][2] = a[3][3] * a[1][1] - a[3][1] * a[1][3];
                            a[3][2] = a[1][3] * a[2][1] - a[1][1] * a[2][3];
                        }//end if
                    }
                }//30
                if (mode != 0) {
                    if (mode != zero) {
                        for (int l = 1; l <= 2; l++) {//30 do l = 1, 2
                            d = zero;
                            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                                b[i][l] = r[i][1] * a[1][l] + r[i][2] * a[2][l] + r[i][3] * a[3][l];
                                d = d + Math.pow(b[i][l], 2);
                            }//end do
                            if (d > zero) d = 1.0 / Math.sqrt(d);
                            for (int i = 1; i <= 3; i++) {//do i = 1, 3
                                b[i][l] = b[i][l] * d;
                            }//end do
                        }//end do
                    }
                    d = b[1][1] * b[1][2] + b[2][1] * b[2][2] + b[3][1] * b[3][2];
                    p = zero;

                    for (int i = 1; i <= 3; i++) {//do i = 1, 3
                        b[i][2] = b[i][2] - d * b[i][1];
                        p = p + Math.pow(b[i][2], 2);
                    }//end do
                    if (p <= tol) {//then
                        p = 1.0;
                        for (int i = 1; i <= 3; i++) {//do i = 1, 3
                            if (p < Math.abs(b[i][1])) continue; //cycle
                            p = Math.abs(b[i][1]);
                            j = i;
                        }//end do
                        k = ip2312[j];
                        l = ip2312[j + 1];
                        p = Math.sqrt(b[k][1] * b[k][1] + b[l][1] * b[l][1]);
                        if (p > tol) { //if (p <= tol) goto 40
                            b[j][2] = zero;
                            b[k][2] = -b[l][1] / p;
                            b[l][2] = b[k][1] / p;
                        }
                    } else {
                        p = 1.0 / Math.sqrt(p);
                        //do i = 1, 3;
                        for (int i = 1; i <= 3; i++) {
                            b[i][2] = b[i][2] * p;
                        }//end do
                    }//end if
                    if (p > tol) {
                        b[1][3] = b[2][1] * b[3][2] - b[2][2] * b[3][1];
                        b[2][3] = b[3][1] * b[1][2] - b[3][2] * b[1][1];
                        b[3][3] = b[1][1] * b[2][2] - b[1][2] * b[2][1];

                        for (int i = 1; i <= 3; i++) {//do i = 1, 3
                            for (int j = 1; j <= 3; j++) {//do j = 1, 3
                                u[i][j] = b[i][1] * a[j][1] + b[i][2] * a[j][2] + b[i][3] * a[j][3];
                            }//end do
                        }//end do
                    }
                } // 40
                //40 do i = 1, 3
                for (int i = 1; i <= 3; i++) {
                    t[i] = ((yc[i] - u[i][1] * xc[1]) - u[i][2] * xc[2]) - u[i][3] * xc[3];
                }//end do
            }
            //50 do i = 1][3
            for (int i = 1; i <= 3; i++) {
                if (e[i] < zero) e[i] = zero;
                e[i] = Math.sqrt(e[i]);
            }//end do

            ier = 0;
            if (e[2] <= (e[1] * 1.0d - 05)) ier = -1;
            d = e[3];
            if (sigma < 0.0) {//then
                d = -d;
                if ((e[2] - e[3]) <= (e[1] * 1.0d - 05)) ier = -1;
            }//end if
            d = (d + e[2]) + e[1];

            rms = (e0 - d) - d;
            if (rms < 0.0) rms = 0.0;

            return;
            // end
        }
    }

    /*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    c     1, collect those residues with dis<d;
    c     2, calculate score_GDT, score_maxsub, score_TM
    ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   */
    private static class ScoreFun {//subroutine score_fun
        static final int nmax = 3000; //PARAMETER(nmax=3000)

        //common/stru/

        //double[] xa = new double[nmax], ya = new double[nmax];
        //double[] za = new double[nmax], xb = new double[nmax], yb = new double[nmax], zb = new double[nmax];
        //common/nres/

        //int[] nresA = new int[nmax], nresB = new int[nmax];
        //int nseqA, nseqB;
        double d;//common/para/d,d0
                    /*common/align/n_ali, */

        //int[] iA = new int[nmax], iB = new int[nmax];
        //common/nscore/

        //int[] i_ali = new int[nmax];

        //common/scores/score
        double score;

        double d_tmp = d;

        //21
        public ScoreFun(boolean debug) {
            double score_sum = 0;
            n_cut = 0; //number of residue-pairs dis<d, for iteration
            while ((n_cut <= 3) & (n_ali > 3)) {
                n_cut = 0; //number of residue-pairs dis<d, for iteration
                score_sum = 0; //               !TMscore
                for (int k = 1; k <= n_ali; k++) {//do k=1,n_ali
                    int i = iA[k];                //![1,nseqA] reoder number of structureA
                    int j = iB[k]; //               ![1,nseqB]
                    if ((i != 0) & (j != 0)) {
                        double dx = xa[i] - xb[j];
                        double dy = ya[i] - yb[j];
                        double dz = za[i] - zb[j];
                        double dis = Math.sqrt(dx * dx + dy * dy + dz * dz);
                        if (dis < d_tmp) { //then
                            n_cut = n_cut + 1;
                            i_ali[n_cut] = k;//      ![1,n_ali], mark the residue-pairs in dis<d
                        }//endif
                        //Utils.printDebug(this,"cccccccccc "+d0+" "+i+" "+j+" "+dis+" "+(dis / d0)+" "+(1 / (1 + (dis / d0) * (dis / d0))));
                        score_sum = score_sum + 1 / (1 + (dis / d0) * (dis / d0));
                    }
                }//enddo
                if ((n_cut <= 3) & (n_ali > 3)) { //then
                    d_tmp = d_tmp + 0.5;
                }
                //goto 21
            }//endif
            score = score_sum / nseqB; //!TM-score

            return;
            // end
        }
    }
}
                //--------------------------------------------------------

/*  http://zhanglab.ccmb.med.umich.edu/TM-score/TMscore_subroutine.f
* **************************************************************************
        *     This is a simple example for showing how to call the TM-score
        *     subroutine to calculate the TM-score superpositions.
        **************************************************************************

        program example
        PARAMETER(nmax=3000)
        dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
        dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

        ********* read the structure1--------->
        open(unit=10,file='pdb1',status='old')
        read(10,*)L1
        do i=1,L1
        read(10,103)du,n1(i),du,x1(i),y1(i),z1(i)
        enddo
        close(10)
        ********* read the structure2--------->
        open(unit=10,file='pdb2',status='old')
        read(10,*)L2
        do i=1,L2
        read(10,103)du,n2(i),du,x2(i),y2(i),z2(i)
        enddo
        close(10)
        103  format(A22,i4,A4,3F8.3)

        ********* calculate TM-score ---------->
        call TMscore(L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,TM,Rcomm,Lcomm)
        do i=1,L1
        write(*,*)i,n1(i),x1(i),y1(i),z1(i)
        enddo
        do i=1,L2
        write(*,*)i,n2(i),x2(i),y2(i),z2(i)
        enddo
        write(*,*)'TMscore=',TM
        write(*,*)'Number of residues in common=',Lcomm
        write(*,*)'RMSD of the common residues=',Rcomm

        stop
        end

        *************************************************************************
        *************************************************************************
        *     This is a subroutine to compare two structures and find the
        *     superposition that has the maximum TM-score.
        *     Reference: Yang Zhang, Jeffrey Skolnick, Proteins 2004 57:702-10.
        *
        *     Explanations:
        *     L1--Length of the first structure
        *     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
        *     n1(i)--Residue sequence number of i'th residue at the first structure
        *     L2--Length of the second structure
        *     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
        *     n2(i)--Residue sequence number of i'th residue at the second structure
        *     TM--TM-score of the comparison
        *     Rcomm--RMSD of two structures in the common aligned residues
        *     Lcomm--Length of the common aligned regions
        *
        *     Note:
        *     1, Always put native as the second structure, by which TM-score
        *        is normalized.
        *     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
        *        TM-score superposition.
        *************************************************************************
        *************************************************************************
        subroutine TMscore(L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,TM,Rcomm,Lcomm)
        PARAMETER(nmax=3000)
        common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
        common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
        common/para/d,d0
        common/align/n_ali,iA(nmax),iB(nmax)
        common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
        dimension k_ali(nmax),k_ali0(nmax)
        dimension L_ini(100),iq(nmax)
        common/scores/score
        double precision score,score_max
        dimension xa(nmax),ya(nmax),za(nmax)

        dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
        dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

        ccc   RMSD:
        double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
        double precision u(3,3),t(3),rms,drms !armsd is real
        data w /nmax*1.0/
        ccc

        ********* convert input data ****************
        nseqA=L1
        do i=1,nseqA
        xa(i)=x1(i)
        ya(i)=y1(i)
        za(i)=z1(i)
        nresA(i)=n1(i)
        enddo
        nseqB=L2
        do i=1,L2
        xb(i)=x2(i)
        yb(i)=y2(i)
        zb(i)=z2(i)
        nresB(i)=n2(i)
        enddo

        ******************************************************************
        *     pickup the aligned residues:
        ******************************************************************
        k=0
        do i=1,nseqA
           do j=1,nseqB
              if(nresA(i).eq.nresB(j))then
                 k=k+1
                 iA(k)=i
                 iB(k)=j
                goto 205
              endif
           enddo
           205     continue
        enddo
        n_ali=k                   !number of aligned residues
        Lcomm=n_ali
        if(n_ali.lt.1)then
        c         write(*,*)'There is no common residues in the input structures'
        TM=0
        Rcomm=0
        return
        endif

        ************/////
        /*     parameters:
        *****************
        ***   d0------------->
        d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
        if(d0.lt.0.5)d0=0.5
        ***   d0_search ----->
        d0_search=d0
        if(d0_search.gt.8)d0_search=8
        if(d0_search.lt.4.5)d0_search=4.5
        ***   iterative parameters ----->
        n_it=20                   !maximum number of iterations
        d_output=5                !for output alignment
        n_init_max=6              !maximum number of L_init
        n_init=0
        L_ini_min=4
        if(n_ali.lt.4)L_ini_min=n_ali
        do i=1,n_init_max-1
        n_init=n_init+1
        L_ini(n_init)=n_ali/2**(n_init-1)
        if(L_ini(n_init).le.L_ini_min)then
        L_ini(n_init)=L_ini_min
        goto 402
        endif
        enddo
        n_init=n_init+1
        L_ini(n_init)=L_ini_min
        402  continue

        ******************************************************************
        *     find the maximum score starting from local structures superposition
        ******************************************************************
        score_max=-1              !TM-score
        do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
        LL=0
        ka=0
        do i=1,L_init
        k=iL+i-1          ![1,n_ali] common aligned
        r_1(1,i)=xa(iA(k))
        r_1(2,i)=ya(iA(k))
        r_1(3,i)=za(iA(k))
        r_2(1,i)=xb(iB(k))
        r_2(2,i)=yb(iB(k))
        r_2(3,i)=zb(iB(k))
        LL=LL+1
        ka=ka+1
        k_ali(ka)=k
        enddo
        call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
        if(i_init.eq.1)then  !global superposition
        armsd=dsqrt(rms/LL)
        Rcomm=armsd
        endif
        do j=1,nseqA
        xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
        yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
        zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
        enddo
        d=d0_search-1
        call score_fun       !init, get scores, n_cut+i_ali(i) for iteration
        if(score_max.lt.score)then
        score_max=score
        ka0=ka
        do i=1,ka0
        k_ali0(i)=k_ali(i)
        enddo
        endif
        ***   iteration for extending ---------------------------------->
        d=d0_search+1
        do 301 it=1,n_it
        LL=0
        ka=0
        do i=1,n_cut
        m=i_ali(i)     ![1,n_ali]
        r_1(1,i)=xa(iA(m))
        r_1(2,i)=ya(iA(m))
        r_1(3,i)=za(iA(m))
        r_2(1,i)=xb(iB(m))
        r_2(2,i)=yb(iB(m))
        r_2(3,i)=zb(iB(m))
        ka=ka+1
        k_ali(ka)=m
        LL=LL+1
        enddo
        call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
        do j=1,nseqA
        xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
        yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
        zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
        enddo
        call score_fun    !get scores, n_cut+i_ali(i) for iteration
        if(score_max.lt.score)then
        score_max=score
        ka0=ka
        do i=1,ka
        k_ali0(i)=k_ali(i)
        enddo
        endif
        if(it.eq.n_it)goto 302
        if(n_cut.eq.ka)then
        neq=0
        do i=1,n_cut
        if(i_ali(i).eq.k_ali(i))neq=neq+1
        enddo
        if(n_cut.eq.neq)goto 302
        endif
        301       continue             !for iteration
        302       continue
        300    continue                !for shift
        333  continue                  !for initial length, L_ali/M

        ******** return the final rotation ****************
        LL=0
        do i=1,ka0
        m=k_ali0(i)            !record of the best alignment
        r_1(1,i)=xa(iA(m))
        r_1(2,i)=ya(iA(m))
        r_1(3,i)=za(iA(m))
        r_2(1,i)=xb(iB(m))
        r_2(2,i)=yb(iB(m))
        r_2(3,i)=zb(iB(m))
        LL=LL+1
        enddo
        call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
        do j=1,nseqA
        x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
        y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
        z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
        enddo
        TM=score_max

        c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        return
        END

        ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
        ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine score_fun
        PARAMETER(nmax=3000)

        common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
        common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
        common/para/d,d0
        common/align/n_ali,iA(nmax),iB(nmax)
        common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
        common/scores/score
        double precision score,score_max

        d_tmp=d
        21   n_cut=0                   !number of residue-pairs dis<d, for iteration
        score_sum=0               !TMscore
        do k=1,n_ali
        i=iA(k)                ![1,nseqA] reoder number of structureA
        j=iB(k)                ![1,nseqB]
        dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
        if(dis.lt.d_tmp)then
        n_cut=n_cut+1
        i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
                 score_sum=score_sum+1/(1+(dis/d0)**2)
                 enddo
                 if(n_cut.lt.3.and.n_ali.gt.3)then
                 d_tmp=d_tmp+.5
                 goto 21
                 endif
                 score=score_sum/float(nseqB) !TM-score

                 return
                 end

                 cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
                 c  w    - w(m) is weight for atom pair  c m           (given)
                 c  x    - x(i,m) are coordinates of atom c m in set x       (given)
                 c  y    - y(i,m) are coordinates of atom c m in set y       (given)
                 c  n    - n is number of atom pairs                         (given)
                 c  mode  - 0:calculate rms only                             (given)
                 c          1:calculate rms,u,t                              (takes longer)
                 c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
                 c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
                 c  t    - t(i)   is translation vector for best superposition  (result)
                 c  ier  - 0: a unique optimal superposition has been determined(result)
                 c       -1: superposition is not unique but optimal
                 c       -2: no result obtained because of negative weights w
                 c           or all weights equal to zero.
                 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                 subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
                 double precision w(*), x(3,*), y(3,*)
                 integer n, mode

                 double precision rms, u(3,3), t(3)
                 integer ier

                 integer i, j, k, l, m1, m
                 integer ip(9), ip2312(4)
                 double precision r(3,3), xc(3), yc(3), wc
                 double precision a(3,3), b(3,3), e(3), rr(6), ss(6)
                 double precision e0, d, spur, det, cof, h, g
                 double precision cth, sth, sqrth, p, sigma

                 double precision sqrt3, tol, zero

                 data sqrt3 / 1.73205080756888d+00 /
                 data tol / 1.0d-2 /
                 data zero / 0.0d+00 /
                 data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
                 data ip2312 / 2, 3, 1, 2 /

                 wc = zero
                 rms = zero
                 e0 = zero

                 do i=1, 3
                 xc(i) = zero
                 yc(i) = zero
                 t(i) = zero
                 do j=1, 3
                 r(i,j) = zero
                 u(i,j) = zero
                 a(i,j) = zero
                 if( i .eq. j ) then
                 u(i,j) = 1.0
                 a(i,j) = 1.0
                 end if
                 end do
                 end do

                 ier = -1
                 if( n .lt. 1 ) return
                 ier = -2
                 do m=1, n
                 if( w(m) .lt. 0.0 ) return
                 wc = wc + w(m)
                 do i=1, 3
                 xc(i) = xc(i) + w(m)*x(i,m)
                 yc(i) = yc(i) + w(m)*y(i,m)
                 end do
                 end do
                 if( wc .le. zero ) return
                 do i=1, 3
                 xc(i) = xc(i) / wc
                 yc(i) = yc(i) / wc
                 end do

                 do m=1, n
                 do i=1, 3
                 e0=e0+w(m)*((x(i,m)-xc(i))**2+(y(i,m)-yc(i))**2)
                 d = w(m) * ( y(i,m) - yc(i) )
                 do j=1, 3
                 r(i,j) = r(i,j) + d*( x(j,m) - xc(j) )
                 end do
                 end do
                 end do

                 det = r(1,1) * ( (r(2,2)*r(3,3)) - (r(2,3)*r(3,2)) )
                 &     - r(1,2) * ( (r(2,1)*r(3,3)) - (r(2,3)*r(3,1)) )
                 &     + r(1,3) * ( (r(2,1)*r(3,2)) - (r(2,2)*r(3,1)) )

                 sigma = det

                 m = 0;
                 do j=1, 3
                 do i=1, j
                 m = m+1
                 rr(m) = r(1,i)*r(1,j) + r(2,i)*r(2,j) + r(3,i)*r(3,j)
                 end do
                 end do

                 spur = (rr(1)+rr(3)+rr(6)) / 3.0
                 cof = (((((rr(3)*rr(6) - rr(5)*rr(5)) + rr(1)*rr(6))
                 &     - rr(4)*rr(4)) + rr(1)*rr(3)) - rr(2)*rr(2)) / 3.0
                 det = det*det

                 do i=1, 3
                 e(i) = spur
                 end do
                 if( spur .le. zero ) goto 40
                 d = spur*spur
                 h = d - cof
                 g = (spur*cof - det)/2.0 - spur*h
                 if( h .le. zero ) then
                 if( mode .eq. 0 ) then
                 goto 50
                 else
                 goto 30
                 end if
                 end if
                 sqrth = dsqrt(h)
                 d = h*h*h - g*g
                 if( d .lt. zero ) d = zero
                 d = datan2( dsqrt(d), -g ) / 3.0
                 cth = sqrth * dcos(d)
                 sth = sqrth*sqrt3*dsin(d)
                 e(1) = (spur + cth) + cth
                 e(2) = (spur - cth) + sth
                 e(3) = (spur - cth) - sth

                 if( mode .eq. 0 ) then
                 goto 50
                 end if

                 do l=1, 3, 2
                 d = e(l)
                 ss(1) = (d-rr(3)) * (d-rr(6))  - rr(5)*rr(5)
                 ss(2) = (d-rr(6)) * rr(2)      + rr(4)*rr(5)
                 ss(3) = (d-rr(1)) * (d-rr(6))  - rr(4)*rr(4)
                 ss(4) = (d-rr(3)) * rr(4)      + rr(2)*rr(5)
                 ss(5) = (d-rr(1)) * rr(5)      + rr(2)*rr(4)
                 ss(6) = (d-rr(1)) * (d-rr(3))  - rr(2)*rr(2)

                 if( dabs(ss(1)) .ge. dabs(ss(3)) ) then
                 j=1
                 if( dabs(ss(1)) .lt. dabs(ss(6)) ) j = 3
                 else if( dabs(ss(3)) .ge. dabs(ss(6)) ) then
                 j = 2
                 else
                 j = 3
                 end if

                 d = zero
                 j = 3 * (j - 1)

                 do i=1, 3
                 k = ip(i+j)
                 a(i,l) = ss(k)
                 d = d + ss(k)*ss(k)
                 end do
                 if( d .gt. zero ) d = 1.0 / dsqrt(d)
                 do i=1, 3
                 a(i,l) = a(i,l) * d
                 end do
                 end do

                 d = a(1,1)*a(1,3) + a(2,1)*a(2,3) + a(3,1)*a(3,3)
                 if ((e(1) - e(2)) .gt. (e(2) - e(3))) then
                 m1 = 3
                 m = 1
                 else
                 m1 = 1
                 m = 3
                 endif

                 p = zero
                 do i=1, 3
                 a(i,m1) = a(i,m1) - d*a(i,m)
                 p = p + a(i,m1)**2
                 end do
                 if( p .le. tol ) then
                 p = 1.0
                 do i=1, 3
                 if (p .lt. dabs(a(i,m))) cycle
                 p = dabs( a(i,m) )
                 j = i
                 end do
                 k = ip2312(j)
                 l = ip2312(j+1)
                 p = dsqrt( a(k,m)**2 + a(l,m)**2 )
                 if( p .le. tol ) goto 40
                 a(j,m1) = zero
                 a(k,m1) = -a(l,m)/p
                 a(l,m1) =  a(k,m)/p
                 else
                 p = 1.0 / dsqrt(p)
                 do i=1, 3
                 a(i,m1) = a(i,m1)*p
                 end do
                 end if

                 a(1,2) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
                 a(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
                 a(3,2) = a(1,3)*a(2,1) - a(1,1)*a(2,3)

                 30   do l=1, 2
                 d = zero
                 do i=1, 3
                 b(i,l) = r(i,1)*a(1,l) + r(i,2)*a(2,l) + r(i,3)*a(3,l)
                 d = d + b(i,l)**2
                 end do
                 if( d .gt. zero ) d = 1.0 / dsqrt(d)
                 do i=1, 3
                 b(i,l) = b(i,l)*d
                 end do
                 end do
                 d = b(1,1)*b(1,2) + b(2,1)*b(2,2) + b(3,1)*b(3,2)
                 p = zero

                 do i=1, 3
                 b(i,2) = b(i,2) - d*b(i,1)
                 p = p + b(i,2)**2
                 end do
                 if( p .le. tol ) then
                 p = 1.0
                 do i=1, 3
                 if( p .lt. dabs(b(i,1)) ) cycle
                 p = dabs( b(i,1) )
                 j = i
                 end do
                 k = ip2312(j)
                 l = ip2312(j+1)
                 p = dsqrt( b(k,1)**2 + b(l,1)**2 )
                 if( p .le. tol ) goto 40
                 b(j,2) = zero
                 b(k,2) = -b(l,1)/p
                 b(l,2) =  b(k,1)/p
                 else
                 p = 1.0 / dsqrt(p)
                 do i=1, 3
                 b(i,2) = b(i,2)*p
                 end do
                 end if

                 b(1,3) = b(2,1)*b(3,2) - b(2,2)*b(3,1)
                 b(2,3) = b(3,1)*b(1,2) - b(3,2)*b(1,1)
                 b(3,3) = b(1,1)*b(2,2) - b(1,2)*b(2,1)

                 do i=1, 3
                 do j=1, 3
                 u(i,j) = b(i,1)*a(j,1) + b(i,2)*a(j,2) + b(i,3)*a(j,3)
                 end do
                 end do

                 40   do i=1, 3
                 t(i) = ((yc(i) - u(i,1)*xc(1)) - u(i,2)*xc(2)) - u(i,3)*xc(3)
                 end do
                 50   do i=1, 3
                 if( e(i) .lt. zero ) e(i) = zero
                 e(i) = dsqrt( e(i) )
                 end do

                 ier = 0
                 if( e(2) .le. (e(1) * 1.0d-05) ) ier = -1

                 d = e(3)
                 if( sigma .lt. 0.0 ) then
                 d = - d
                 if( (e(2) - e(3)) .le. (e(1) * 1.0d-05) ) ier = -1
                 end if
                 d = (d + e(2)) + e(1)

                 rms = (e0 - d) - d
                 if( rms .lt. 0.0 ) rms = 0.0

                 return
                 end
*/