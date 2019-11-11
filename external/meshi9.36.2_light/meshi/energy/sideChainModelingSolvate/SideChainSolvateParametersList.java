/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.sideChainModelingSolvate;

import java.io.*;
import java.util.*;

import meshi.util.Utils;
import meshi.util.mathTools.*;
import meshi.energy.*;
import meshi.parameters.*;

/**
 * This class is a truncated form of the corresponding class in the "meshi.energy.solvation" package. It was designed
 * solely for accelerating the application SCMOD (concurrent sidechain modeling), and should not be used for other
 * purposes.
 * <p/>
 * Please do not use it!
 */

public class SideChainSolvateParametersList extends ArrayList implements Parameters {

    public final int NTsai = 14; // The number of atom types used in Tsai 99'. Any Hydrogen is type 14.
    public final int[] atomicTypeConverter; // From MESHI atom types to Tsai 99'.
    public final double maxEnd; // The maximal distance where any sigmoid is not zero. This value must be less than the Rmax in the distance matrix. 
    public final Spline1D[] atomTypeSplines;
    public final double[][] Cend;
    public final double[][] Cp1;
    public final double[][] Cp2;
    public final double[][] CvalAtp1;
    public final double[][] CvalAtp2;
    public final double[][] HBend;
    public final double[][] HBp1;
    public final double[][] HBp2;
    public final double[][] HBvalAtp1;
    public final double[][] HBvalAtp2;

    /**
     * The parameter to the constructor is an array of 12 Strings, giving the 12 files required as parameters.
     * See the MeshiPotential class for a detailed list.
     */
    public SideChainSolvateParametersList(String[] parameterFiles) {
        super();
        int maxAtomType = -1;
        int tmp;
        BufferedReader br;
        StringTokenizer stok;
        String line = "";
        // Converting the 190 MESHI atom types to the 14 mentioned in Tsai 99'
        Utils.println("Reading solvation parameter file: " + parameterFiles[0]);
        try {
            // first pass isOn the file - to find the maximal atom type
            br = new BufferedReader(new FileReader(parameterFiles[0]));
            line = br.readLine();
            while (line != null) {
                stok = new StringTokenizer(line);
                tmp = (AtomType.type(stok.nextToken().trim())).ordinal();
                if (tmp > maxAtomType)
                    maxAtomType = tmp;
                line = br.readLine();
            }
            br.close();
            atomicTypeConverter = new int[maxAtomType + 1];
            for (int c = 0; c < atomicTypeConverter.length; c++)
                atomicTypeConverter[c] = -1;
            // second pass isOn the file - reading the new types
            br = new BufferedReader(new FileReader(parameterFiles[0]));
            line = br.readLine();
            while (line != null) {
                stok = new StringTokenizer(line);
                tmp = (AtomType.type(stok.nextToken().trim())).ordinal();
                atomicTypeConverter[tmp] = Integer.valueOf(stok.nextToken().trim()).intValue() - 1;
                line = br.readLine();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }

        // Reading the spline parameters for the 190 atom types.
        atomTypeSplines = new Spline1D[maxAtomType + 1];
        Utils.println("Reading solvation parameter file: " + parameterFiles[1]);
        try {
            br = new BufferedReader(new FileReader(parameterFiles[1]));
            line = br.readLine();
            while (line != null) {
                stok = new StringTokenizer(line);
                tmp = (AtomType.type(stok.nextToken().trim())).ordinal();
                atomTypeSplines[tmp] = new Spline1D(stok);
                line = br.readLine();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }

        Cend = new double[NTsai][NTsai];
        readSigmoidValueFile(Cend, parameterFiles[2]);

        Cp1 = new double[NTsai][NTsai];
        readSigmoidValueFile(Cp1, parameterFiles[3]);

        Cp2 = new double[NTsai][NTsai];
        readSigmoidValueFile(Cp2, parameterFiles[4]);

        CvalAtp1 = new double[NTsai][NTsai];
        readSigmoidValueFile(CvalAtp1, parameterFiles[5]);

        CvalAtp2 = new double[NTsai][NTsai];
        readSigmoidValueFile(CvalAtp2, parameterFiles[6]);

        HBend = new double[NTsai][NTsai];
        readSigmoidValueFile(HBend, parameterFiles[7]);

        HBp1 = new double[NTsai][NTsai];
        readSigmoidValueFile(HBp1, parameterFiles[8]);

        HBp2 = new double[NTsai][NTsai];
        readSigmoidValueFile(HBp2, parameterFiles[9]);

        HBvalAtp1 = new double[NTsai][NTsai];
        readSigmoidValueFile(HBvalAtp1, parameterFiles[10]);

        HBvalAtp2 = new double[NTsai][NTsai];
        readSigmoidValueFile(HBvalAtp2, parameterFiles[11]);

        // Finding the largest ending distance
        double maxEndTmp = -1.0;
        for (int c1 = 0; c1 < NTsai; c1++)
            for (int c2 = 0; c2 < NTsai; c2++) {
                if (Cend[c1][c2] > maxEndTmp)
                    maxEndTmp = Cend[c1][c2];
                if (HBend[c1][c2] > maxEndTmp)
                    maxEndTmp = HBend[c1][c2];
            }
        maxEnd = maxEndTmp;

    }

    public Parameters createParameters(String s) {
        throw new RuntimeException("This method should not be called");
    }

    public Parameters parameters(Object obj) {
        throw new RuntimeException("This method should not be called");
    }


    // Reading a 14x14 file into the relevent parameter.

    private void readSigmoidValueFile(double[][] ar, String filename) {
        BufferedReader br;
        StringTokenizer stok;
        String line = "";
        Utils.println("Reading solvation parameter file: " + filename);
        try {
            br = new BufferedReader(new FileReader(filename));
            for (int ind1 = 0; ind1 < NTsai; ind1++) {
                line = br.readLine();
                if (line == null)
                    throw new RuntimeException("In " + filename + " there should be a " + NTsai + "x" +
                            NTsai + " array of doubles");
                stok = new StringTokenizer(line);
                for (int ind2 = 0; ind2 < NTsai; ind2++)
                    ar[ind1][ind2] = Double.valueOf(stok.nextToken().trim()).doubleValue();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }
}
