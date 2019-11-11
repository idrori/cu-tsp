/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

/**
 * Coefficients matrix for a single atom pair spline.
 */
public class CoefficientsMatrixForAtomPairSpline {
    public double[] bins;        /* list of bins; later used in calculation */
    public double[][] coefs;
    public final String name;    /* names of atoms these coefficients calculate*/

    /**
     * Creates new coefficients matrix. We assume each break has four coefficients.
     */
    public CoefficientsMatrixForAtomPairSpline(double[] bins, int numBins, String name) {
        this.bins = bins;
        coefs = new double[numBins][4];
        this.name = name;
    }

    /**
     * Calculate value of polynomial spline and first derivative for given distance value.
     */
    public double[] evaluate(double x) {
        double[] energy = new double[2];

        /* find bin */
        int bin = 0;
        while (bin < bins.length - 1 && bins[bin + 1] <= x)
            bin++;

        /* correct x according to its distance from start of bin */
        x = x - bins[bin];
        /* calculate powers of x */
        double x2 = x * x;
        double x3 = x2 * x;

        /* calculate energy */
        energy[0] = coefs[bin][0] * x3 + coefs[bin][1] * x2 + coefs[bin][2] * x + coefs[bin][3];
        /* calculate first derivative */
        energy[1] = 3 * coefs[bin][0] * x2 + 2 * coefs[bin][1] * x + coefs[bin][2];

        return energy;
    }

    /**
     * Print all coefficients to screen.
     */
    public void printAllCoefficients() {
        for (int i = 0; i < coefs.length; i++)
            System.out.println(coefs[i][0] + ", " + coefs[i][1] + ", " + coefs[i][2] + ", " + coefs[i][3] );
	}

}
