/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions;

import meshi.energy.simpleEnergyTerms.*;


/**
 * Polynomial for calculating energy values for a given amino acid, and given
 * torsion angles.
 *
 * @author El-ad David Amir
 */
public class SplinedPolynomial
        implements CompositeTorsionsDefinitions {

    private int numVariables;        /* how many variables the polynom accepts */
    private int aac;            /* amino acid type */
    private int torsionAngles[];    /* torsion angles used */
    private int ssType;                /* secondary structure type */
    private double breaks[][];        /* polynomial breaks for each variable */
    private double coefs[];             /* polynomial coefficients */
    private double[][] quickPowerMatrix = new double[3][4];

    private int order[];        /* order of polynomial */

    /**
     * Creates a new splined polynomial.
     */
    public SplinedPolynomial(
            int numVariables, int aac, int torsionAngles[], int ssType,
            double breaks[][], double coefs[]) {
        this.numVariables = numVariables;
        this.aac = aac;
        this.torsionAngles = torsionAngles;
        this.ssType = ssType;
        this.breaks = breaks;
        this.coefs = coefs;

        /* polynomial's order is assumed to be four */
        order = new int[numVariables];
        for (int i = 0; i < numVariables; i++)
            order[i] = 4;
    }

    /**
     * checks whether this polynomial is of the parameters needed.
     */
    public boolean isPolynomialNeeded(int aac, int torsionAngles[], int ssType) {
        if (this.aac != aac)
            return false;

        if (this.torsionAngles.length != torsionAngles.length)
            return false;

        for (int i = 0; i < torsionAngles.length; i++)
            if (this.torsionAngles[i] != torsionAngles[i])
                return false;

        if (this.ssType != ssType)
            return false;

        return true;
    }

    /**
     * Calculates polynomial's derivation for given variable with
     * value for list of arguments.
     *
     * @param derivVar variable to be derived (zero for calculation of polynomial).
     */
    public double value(int derivVar, double... args) {
        double result = 0;        /* result of calculation */
        int i, j, k;

        /* verify number of arguments is correct */
        if (numVariables != args.length)
            throw new RuntimeException("polynomial number of variables mismatch");
        /* verify we're not asking to derive a non existent variable */
        if (derivVar > numVariables)
            throw new RuntimeException("derived variable identifier doesn't exist in this polynomial");

        /* calculate 1D polynomial */
        if (args.length == 1) {
            double tor1 = args[0];
            int bin1 = findBin(breaks[0], tor1);

            /* verify a bin was found */
            if (bin1 == -1)
                throw new RuntimeException(
                        "unable to locate correct bin for first torsion angle (" +
                                tor1 + "), " + toString());

            tor1 = fixTorsionToBin(breaks[0], tor1, bin1);

            /* convert torsion angle to distance from beginning of bin */
            tor1 = tor1 - breaks[0][bin1];
            quickPowerMatrix[0][0] = 1;
            quickPowerMatrix[0][1] = tor1;
            quickPowerMatrix[0][2] = tor1 * tor1;
            quickPowerMatrix[0][3] = tor1 * tor1 * tor1;
            /* calculation acceleration */
            int bin1x4 = bin1 * 4;

            /* calculation for univariate splined polynomial */
            if (derivVar == 0)
                /* no derivation required */
                for (i = 0; i < 4; i++)
                    result = result +
                            coefs[bin1x4 + i] * quickPowerMatrix[0][4 - (i + 1)];
            else if (derivVar == 1)
                /* derive first variable */
                for (i = 0; i < 3; i++)
                    result = result +
                            coefs[bin1x4 + i] * quickPowerMatrix[0][3 - (i + 1)] * (4 - (i + 1));
        }
        /* calculate 2D polynomial */
        else if (args.length == 2) {
            double tor1 = args[0];
            double tor2 = args[1];
            int bin1 = findBin(breaks[0], tor1);
            int bin2 = findBin(breaks[1], tor2);

            /* verify a bin was found */
            if (bin1 == -1)
                throw new RuntimeException(
                        "unable to locate correct bin for first torsion angle (" +
                                tor1 + "), " + toString());
            if (bin2 == -1)
                throw new RuntimeException(
                        "unable to locate correct bin for second torsion angle (" +
                                tor2 + "), " + toString());

            tor1 = fixTorsionToBin(breaks[0], tor1, bin1);
            tor2 = fixTorsionToBin(breaks[1], tor2, bin2);

            /* convert torsion angles to distance from beginning of bin */
            tor1 = tor1 - breaks[0][bin1];
            tor2 = tor2 - breaks[1][bin2];
            quickPowerMatrix[0][0] = 1;
            quickPowerMatrix[0][1] = tor1;
            quickPowerMatrix[0][2] = tor1 * tor1;
            quickPowerMatrix[0][3] = tor1 * tor1 * tor1;
            quickPowerMatrix[1][0] = 1;
            quickPowerMatrix[1][1] = tor2;
            quickPowerMatrix[1][2] = tor2 * tor2;
            quickPowerMatrix[1][3] = tor2 * tor2 * tor2;

            /* calculation acceleration */
            int bin2x4 = bin2 * 4;
            int co2lenx4 = (breaks[1].length - 1) * 4;
            int bin1xco2lenx4x4 = bin1 * co2lenx4 * 4;

            if (derivVar == 0)
                /* no derivation required */
                for (i = 0; i < 4; i++)
                    for (j = 0; j < 4; j++)
                        result = result +
                                coefs[bin1xco2lenx4x4 + co2lenx4 * i + bin2x4 + j] *
                                        quickPowerMatrix[0][4 - (i + 1)] *
                                        quickPowerMatrix[1][4 - (j + 1)];
            else if (derivVar == 1)
                /* derive according to first variable */
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 4; j++)
                        result = result +
                                coefs[bin1xco2lenx4x4 + co2lenx4 * i + bin2x4 + j] *
                                        quickPowerMatrix[0][3 - (i + 1)] * (4 - (i + 1)) *
                                        quickPowerMatrix[1][4 - (j + 1)];
            else if (derivVar == 2)
                /* derive according to second variable */
                for (i = 0; i < 4; i++)
                    for (j = 0; j < 3; j++)
                        result = result +
                                coefs[bin1xco2lenx4x4 + co2lenx4 * i + bin2x4 + j] *
                                        quickPowerMatrix[0][4 - (i + 1)] *
                                        quickPowerMatrix[1][3 - (j + 1)] * (4 - (j + 1));
        }
        /* calculate 3D polynomial */
        else if (args.length == 3) {
            double tor1 = args[0];
            double tor2 = args[1];
            double tor3 = args[2];
            int bin1 = findBin(breaks[0], tor1);
            int bin2 = findBin(breaks[1], tor2);
            int bin3 = findBin(breaks[2], tor3);

            /* verify a bin was found */
            if (bin1 == -1)
                throw new RuntimeException(
                        "unable to locate correct bin for first torsion angle (" +
                                tor1 + "), " + toString());
            if (bin2 == -1)
                throw new RuntimeException(
                        "unable to locate correct bin for second torsion angle (" +
                                tor2 + "), " + toString());
            if (bin2 == -1)
                throw new RuntimeException(
                        "unable to locate correct bin for third torsion angle (" +
                                tor3 + "), " + toString());

            tor1 = fixTorsionToBin(breaks[0], tor1, bin1);
            tor2 = fixTorsionToBin(breaks[1], tor2, bin2);
            tor3 = fixTorsionToBin(breaks[2], tor3, bin3);

            /* convert torsion angles to distance from beginning of bin */
            tor1 = tor1 - breaks[0][bin1];
            tor2 = tor2 - breaks[1][bin2];
            tor3 = tor3 - breaks[2][bin3];
            quickPowerMatrix[0][0] = quickPowerMatrix[1][0] = quickPowerMatrix[2][0] = 1;
            quickPowerMatrix[0][1] = tor1;
            quickPowerMatrix[0][2] = tor1 * tor1;
            quickPowerMatrix[0][3] = tor1 * tor1 * tor1;
            quickPowerMatrix[1][1] = tor2;
            quickPowerMatrix[1][2] = tor2 * tor2;
            quickPowerMatrix[1][3] = tor2 * tor2 * tor2;
            quickPowerMatrix[2][1] = tor3;
            quickPowerMatrix[2][2] = tor3 * tor3;
            quickPowerMatrix[2][3] = tor3 * tor3 * tor3;

            /* calculation acceleration */
            int bin3x4 = bin3 * 4;
            int co3lenx4 = (breaks[2].length - 1) * 4;
            int bin2xco3lenx16 = bin2 * co3lenx4 * 4;
            int co2lenxco3lenx16 = (breaks[1].length - 1) * co3lenx4 * 4;
            int bin1xco2lenxco3lenx64 = bin1 * co2lenxco3lenx16 * 4;

            /* calculation for tri-variate splined polynomial */
            if (derivVar == 0)
                /* no derivation required */
                for (i = 0; i < 4; i++)
                    for (j = 0; j < 4; j++)
                        for (k = 0; k < 4; k++)
                            result = result +
                                    coefs[bin1xco2lenxco3lenx64 + co2lenxco3lenx16 * i + bin2xco3lenx16 + co3lenx4 * j + bin3x4 + k] *
                                            quickPowerMatrix[0][4 - (i + 1)] *
                                            quickPowerMatrix[1][4 - (j + 1)] *
                                            quickPowerMatrix[2][4 - (k + 1)];
            else if (derivVar == 1)
                /* derive according to first variable */
                for (i = 0; i < 3; i++) // used to be for( i=0; i<4; i++ )
                    for (j = 0; j < 4; j++)
                        for (k = 0; k < 4; k++)
                            result = result +
                                    coefs[bin1xco2lenxco3lenx64 + co2lenxco3lenx16 * i + bin2xco3lenx16 + co3lenx4 * j + bin3x4 + k] *
                                            quickPowerMatrix[0][3 - (i + 1)] * (4 - (i + 1)) *
                                            quickPowerMatrix[1][4 - (j + 1)] *
                                            quickPowerMatrix[2][4 - (k + 1)];
            else if (derivVar == 2)
                /* derive according to second variable */
                for (i = 0; i < 4; i++)
                    for (j = 0; j < 3; j++)
                        for (k = 0; k < 4; k++)
                            result = result +
                                    coefs[bin1xco2lenxco3lenx64 + co2lenxco3lenx16 * i + bin2xco3lenx16 + co3lenx4 * j + bin3x4 + k] *
                                            quickPowerMatrix[0][4 - (i + 1)] *
                                            quickPowerMatrix[1][3 - (j + 1)] * (4 - (j + 1)) *
                                            quickPowerMatrix[2][4 - (k + 1)];
            else if (derivVar == 3)
                /* derive according to third variable */
                for (i = 0; i < 4; i++)
                    for (j = 0; j < 4; j++)
                        for (k = 0; k < 3; k++)
                            result = result +
                                    coefs[bin1xco2lenxco3lenx64 + co2lenxco3lenx16 * i + bin2xco3lenx16 + co3lenx4 * j + bin3x4 + k] *
                                            quickPowerMatrix[0][4 - (i + 1)] *
                                            quickPowerMatrix[1][4 - (j + 1)] *
                                            quickPowerMatrix[2][3 - (k + 1)] * (4 - (k + 1));
        }

        return result;
    }

    /**
     * Find break index (bin) for torsion angle.
     *
     * @return left bound of break, -1 if no break found.
     */
    public int findBin(double breaks[], double torsion) {
        /* update torsion to be in bin range */
        if (torsion < breaks[0])
            torsion = torsion + 2 * Math.PI;
        else if (torsion > breaks[breaks.length - 1])
            torsion = torsion - 2 * Math.PI;

        /* search for torsion in all the bins */
        for (int binCount = 0; binCount < breaks.length; binCount++)
            if (breaks[binCount] <= torsion && breaks[binCount + 1] > torsion)
                return binCount;
        return -1;
    }

    /**
     * Fixes a torsion to be inside its bin. Since some of the breaks
     * are not from -pi to pi (due to the varied grid construction
     * method) we need this addition to verify the polynomial is
     * calculated correctly.
     */
    public double fixTorsionToBin(double breaks[], double torsion, int bin) {
        if (torsion < breaks[bin])
            torsion = torsion + 2 * Math.PI;
        else if (torsion > breaks[bin + 1])
            torsion = torsion - 2 * Math.PI;
        return torsion;
    }

    /**
     * Converts attributes of polynomial to string.
     */
    public String toString() {
        /* convert torsions to string */
        String torsionAnglesStr = "";
        for (int i : torsionAngles)
            torsionAnglesStr = torsionAnglesStr + " " + i;

        /* convert length of breaks to string */
        String breaksStr = "";
        for (double d[] : breaks)
            breaksStr = breaksStr + " " + d.length;

        /* return all attributes */
        return "numVariables: " + numVariables + ", aac: " + aac +
                ", torsionAngles:" + torsionAnglesStr + ", ssType: " + ssType +
                ", breaks:" + breaksStr + ", coefs: " + coefs.length;
    }

    /**
     * Raises a double by the power, with 0 <= power <= 3. For the special
     * exception of power=-1, returns 0 (as it eases polynomial derivation
     * calculations).
     */
    private static double quickPower(double torsion, int power) {
        switch (power) {
            case -1:
                return 0;
            case 0:
                return 1;
            case 1:
                return torsion;
            case 2:
                return torsion * torsion;
            case 3:
                return torsion * torsion * torsion;
        }

        return -1;
    }
}
