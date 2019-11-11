/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions;

import meshi.energy.simpleEnergyTerms.*;

import meshi.util.file.*;

import java.io.*;
import java.util.Vector;

/**
 * Manages the list of splined polynomials. Loads all
 * splined polynomials from parameters data file.
 *
 * @author El-ad David Amir
 */
public class SplinedPolynomialsLoader {

    /* a small padding factor for lefternmost and righternmost break. since our
      * input is text files, there's a small chance the breaks will not form
      * a complete 2*pi circle.
      */
    private static final double BREAKS_PADDING = 1e-8;

    /* vector of all splined polynomials loaded from file */
    private Vector<SplinedPolynomial> splinedPolynomials;

    /**
     * Loads polynomials from given file.
     */
    public SplinedPolynomialsLoader(String fileName) throws IOException {
        /* reads all lines from this file */
        MeshiLineReader mlr = new MeshiLineReader(fileName);

        /* create the vector */
        splinedPolynomials = new Vector<SplinedPolynomial>();

        /* polynomial parameters */
        int numVariables = -1;
        int aac = -1;
        int torsionAngles[];
        int ssType = -1;
        double breaks[][];
        double coefs[];

        int i, j;

        /* start reading polynomials */
        while (mlr.ready()) {
            /* read number of variables in polynomial */
            numVariables = new Integer(mlr.readLine("#")).intValue();

            /* read amino acid type. we decrease one because Matlab
                * enumeration starts with 1, while Java starts with 0. */
            aac = new Integer(mlr.readLine("#")).intValue() - 1;

            /* setResidue up torsion angles array and read torsion angles */
            torsionAngles = new int[numVariables];
            for (i = 0; i < numVariables; i++)
                torsionAngles[i] = new Integer(mlr.readLine()).intValue();

            /* read secondary structure types */
            ssType = new Integer(mlr.readLine()).intValue();

            /* read all breaks arrays */
            breaks = new double[numVariables][];
            int breakSize;
            int numCoefs = 1;
            for (i = 0; i < numVariables; i++) {
                /* read breaks size */
                breakSize = new Integer(mlr.readLine()).intValue();

                /* setResidue up breaks array read contents of break */
                breaks[i] = new double[breakSize];
                for (j = 0; j < breakSize; j++)
                    breaks[i][j] = new Double(mlr.readLine()).doubleValue();

                /* apply padding */
                breaks[i][0] = breaks[i][0] - BREAKS_PADDING;
                breaks[i][breaks[i].length - 1] =
                        breaks[i][breaks[i].length - 1] + BREAKS_PADDING;

                /* update number of coefficients */
                numCoefs = numCoefs * (breakSize - 1) * 4;
            }

            /* verify breaks are correctly aligned */
            for (i = 0; i < numVariables; i++) {
                if (breaks[i][breaks[i].length - 1] - breaks[i][0] < 2 * Math.PI) {
                    /* report incorrect polynomial information */
                    String torsionsStr = "";
                    for (j = 0; j < numVariables; j++)
                        torsionsStr = torsionsStr + torsionAngles[j] + " ";
                    String breaksStr = "";
                    for (j = 0; j < breaks[i].length; j++)
                        breaksStr = breaksStr + breaks[i][j] + " ";
                    throw new RuntimeException(
                            "polynomial breaks do not match for amino acid " +
                                    aac + ", torsion angles are " + torsionsStr +
                                    ", breaks are " + breaksStr);
                }
            }


            /* setResidue up coefficients array and read all coefficients */
            coefs = new double[numCoefs];
            for (i = 0; i < numCoefs; i++)
                coefs[i] = new Double(mlr.readLine()).doubleValue();

            /* add polynomial to vector */
            splinedPolynomials.add(new SplinedPolynomial(
                    numVariables, aac, torsionAngles, ssType, breaks, coefs));
        }
    }

    /**
     * Finds polynomial of given type in vector.
     */
    public SplinedPolynomial findPolynomial(
            int aac, int torsionAngles[], int ssType) {

        /* scan all polynomials in database for this polynomial */
        for (SplinedPolynomial splinedPolynomial : splinedPolynomials)
            if (splinedPolynomial.isPolynomialNeeded(aac, torsionAngles, ssType))
                return splinedPolynomial;

        return null;
    }

    /**
     * report list of polynomials
     */
    public void reportPolynomialsList() {
        for (SplinedPolynomial splinedPolynomial : splinedPolynomials)
            System.out.println(splinedPolynomial );
	}
}


