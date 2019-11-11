/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.parameters.*;
import meshi.util.Utils;
import meshi.util.mathTools.Spline1D;
import meshi.util.mathTools.Spline2D;

/**
 * The parameter class required for the SolvateEnergy class. This class reads the parameters from 14 different files.
 * The 14 file names are given as a String array to the constructor.
 * <p/>
 * The first 10 files, are files that read parameters for the different sigmoid functions used in the calculations
 * of the carbon and HB indices in SolvateEnergy. See the documentation of the Sigma class to see what are the
 * parameters required to define a sigmoid. The formats of all these files are similar. Each file contains a
 * 14x14 matrix of doubles (14 is the number of types in Tsai 99'. They are: 1-PRO-N ; 2-backbone-N ; 3-ASN-ND ; 4-LYS-NZ ;
 * 5-backbone-O ; 6-ASP-OD ; 7-PHE-CG ; 8-PHE-CD ; 9-VAL-CB ; 10-backbone-CA ; 11-ALA-CB ; 12-MET-SD ; 13-CYS-SG ; 14-hydrogens).
 * Value i,j in the matrix, is the relevent sigmoid property of atom type j isOn atom type i. For example, the sixth value in the
 * second row, corresponds to a sigmoid of hydroxyl oxygen (type 6) isOn a backbone nitrogen (type 2).
 * <p/>
 * The 13 file types:
 * ------------------
 * 1) SolvateExtResCend.dat - A 14x14 matrix of doubles. These values are the 'end' values in the sigmoid that
 * calculates the carbon index.
 * <p/>
 * 2) SolvateExtResCp1.dat - A 14x14 matrix of doubles. These values are the 'p1' values in the sigmoid that
 * calculates the carbon index.
 * <p/>
 * 3) SolvateExtResCp2.dat - A 14x14 matrix of doubles. These values are the 'p2' values in the sigmoid that
 * calculates the carbon index.
 * <p/>
 * 4) SolvateExtResCvalAtp1.dat - A 14x14 matrix of doubles. These values are the 'valAtp1' values in the sigmoid that
 * calculates the carbon index.
 * <p/>
 * 5) SolvateExtResCvalAtp2.dat - A 14x14 matrix of doubles. These values are the 'valAtp2' values in the sigmoid that
 * calculates the carbon index.
 * <p/>
 * 6) SolvateExtResHBend.dat - A 14x14 matrix of doubles. These values are the 'end' values in the sigmoid that
 * calculates the HB index.
 * <p/>
 * 7) SolvateExtResHBp1.dat - A 14x14 matrix of doubles. These values are the 'p1' values in the sigmoid that
 * calculates the HB index.
 * <p/>
 * 8) SolvateExtResHBp2.dat - A 14x14 matrix of doubles. These values are the 'p2' values in the sigmoid that
 * calculates the HB index.
 * <p/>
 * 9) SolvateExtResHBvalAtp1.dat - A 14x14 matrix of doubles. These values are the 'valAtp1' values in the sigmoid that
 * calculates the HB index.
 * <p/>
 * 10) SolvateExtResHBvalAtp2.dat - A 14x14 matrix of doubles. These values are the 'valAtp2' values in the sigmoid that
 * calculates the HB index.
 * <p/>
 * 11) SolvateMESHI2Tsai.dat - In MESHI there are 190 defined atom types. In the solvation energy, we sometime use
 * a reduce representation of 14 atom types as defined in Tsai et. al. 99'. The format of this file is:
 * {MESHI atom type (string)} {Tsai type number (int)}
 * {MESHI atom type (string)} {Tsai type number (int)}
 * {MESHI atom type (string)} {Tsai type number (int)}
 * ...
 * {MESHI atom type (string)} {Tsai type number (int)}
 * <p/>
 * 12) SolvateSCpolarSplines.dat - 2D spline parameters for the solvation of side-chain polar getAtoms.
 * A spline is created for each of the of the 190 atom types, but it is a zero spline for those that are
 * not side-chain polar getAtoms. The format of this file is:
 * TRN {2D spline initialization string for TRN (see Spline2D.java)}
 * TRC {zero spline}
 * TRO {2D spline initialization string for TRO (see Spline2D.java)}
 * AH {zero spline)
 * AN {zero spline)
 * ...
 * CSG {2D spline initialization string for CSG (see Spline2D.java)}
 * ...
 * YOH {2D spline initialization string for YOH (see Spline2D.java)}
 * <p/>
 * 13) SolvateSCcarbonSplines.dat - 1D spline parameters for the solvation of side-chain NON-polar getAtoms.
 * A spline is created for each of the of the 190 atom types, but it is a zero spline for those that are
 * not side-chain NON-polar getAtoms. The format of this file is:
 * TRN {zero spline}
 * TRC {zero spline}
 * TRO {zero spline}
 * AH {zero spline}
 * ...
 * DCG {1D spline initialization string for DCG (see Spline1D.java)}
 * ...
 * YCZ {1D spline initialization string for YCZ (see Spline1D.java)}
 * YOH {zero spline}
 * <p/>
 * 14) SolvateBBpolarSplines.dat - 2D spline parameters for the solvation of backbone polar getAtoms.
 * A spline is created for each of the of the 190 atom types, but it is a zero spline for those that are
 * not backbone polar getAtoms. The format of this file is:
 * TRN {zero spline}
 * TRC {zero spline}
 * TRO {zero spline}
 * AH {zero spline)
 * AN {2D spline initialization string for AN (see Spline2D.java)}
 * ...
 * AO {2D spline initialization string for AO (see Spline2D.java)}
 * ...
 * YOH {zero spline}
 */

public class SolvateParametersList implements Parameters {

    public final static int N_TSAI = 14; // The number of atom types used in Tsai 99'. Any Hydrogen is type 14.
    public final int[] atomicTypeConverter; // From MESHI atom types to Tsai 99'.
    public final double maxEnd; // The maximal distance where any sigmoid is not zero. This value must be less than the Rmax in the distance matrix. 
    public final Spline2D[] splines;
    public final double[][] cEnd;
    public final double[][] cP1;
    public final double[][] cP2;
    public final double[][] cValAtp1;
    public final double[][] cValAtp2;
    public final double[][] hbEnd;
    public final double[][] hbP1;
    public final double[][] hbP2;
    public final double[][] hbValAtp1;
    public final double[][] hbValAtp2;
    public static final SigmaParameters[][] sigmaParameters = new SigmaParameters[N_TSAI][N_TSAI];
    public static final SigmaParameters[][] hbParameters = new SigmaParameters[N_TSAI][N_TSAI];

    /**
     * The parameter to the constructor is an array of 14 Strings, giving the 14 files required as parameters.
     * See the MeshiPotential class for a detailed list.
     */
    public SolvateParametersList(String[] parameterFiles) {
        super();
        int maxAtomType = -1;
        int atomTypeOrdinal;
        BufferedReader br;
        StringTokenizer stok;
        String line = "";

        cEnd = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(cEnd, parameterFiles[0]);

        cP1 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(cP1, parameterFiles[1]);

        cP2 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(cP2, parameterFiles[2]);

        cValAtp1 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(cValAtp1, parameterFiles[3]);

        cValAtp2 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(cValAtp2, parameterFiles[4]);

        hbEnd = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(hbEnd, parameterFiles[5]);

        hbP1 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(hbP1, parameterFiles[6]);

        hbP2 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(hbP2, parameterFiles[7]);

        hbValAtp1 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(hbValAtp1, parameterFiles[8]);

        hbValAtp2 = new double[N_TSAI][N_TSAI];
        readSigmoidValueFile(hbValAtp2, parameterFiles[9]);

        // Finding the largest ending distance
        double maxEndTmp = -1.0;
        for (int c1 = 0; c1 < N_TSAI; c1++)
            for (int c2 = 0; c2 < N_TSAI; c2++) {
                if (cEnd[c1][c2] > maxEndTmp)
                    maxEndTmp = cEnd[c1][c2];
                if (hbEnd[c1][c2] > maxEndTmp)
                    maxEndTmp = hbEnd[c1][c2];
            }
        maxEnd = maxEndTmp;

        // Converting the 190 MESHI atom types to the 14 mentioned in Tsai 99'
        Utils.println("Reading solvation parameter file: " + parameterFiles[10]);
        try {
            // first pass isOn the file - to find the maximal atom type
            br = new BufferedReader(new FileReader(parameterFiles[10]));
            line = br.readLine();
            while (line != null) {
                stok = new StringTokenizer(line);
                atomTypeOrdinal = AtomType.type(stok.nextToken().trim()).ordinal();
                if (atomTypeOrdinal > maxAtomType)
                    maxAtomType = atomTypeOrdinal;
                line = br.readLine();
            }
            br.close();
            atomicTypeConverter = new int[maxAtomType + 1];
            for (int c = 0; c < atomicTypeConverter.length; c++)
                atomicTypeConverter[c] = -1;
            // second pass isOn the file - reading the new types
            br = new BufferedReader(new FileReader(parameterFiles[10]));
            line = br.readLine();
            while (line != null) {
                stok = new StringTokenizer(line);
                atomTypeOrdinal = AtomType.type(stok.nextToken().trim()).ordinal();
                atomicTypeConverter[atomTypeOrdinal] = Integer.valueOf(stok.nextToken().trim()).intValue() - 1;
                line = br.readLine();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }

        // Reading the spline parameters for the solvation of side-chain polar groups (There should
        // be one spline for each residue type).
        splines = new Spline2D[maxAtomType + 1];
        Utils.println("Reading solvation parameter file: " + parameterFiles[11]);
        try {
            br = new BufferedReader(new FileReader(parameterFiles[11]));
            line = br.readLine();
            while (line != null) {
                stok = new StringTokenizer(line);
                atomTypeOrdinal = AtomType.type(stok.nextToken().trim()).ordinal();
                splines[atomTypeOrdinal] = new Spline2D(stok);
                line = br.readLine();
            }
            br.close();
        }
        catch (Exception e) {
            System.out.println("********** error ********");
            System.out.println("e.getClass() = "+e.getClass()+"    e.getMessage() = "+e.getMessage());
            System.out.println("*********** error end ************");
            e.printStackTrace();
            throw new RuntimeException("\nThis is weird.");
        }


        for (int tsaiType1 = 0; tsaiType1 < N_TSAI; tsaiType1++)
            for (int tsaiType2 = 0; tsaiType2 < N_TSAI; tsaiType2++) {
                sigmaParameters[tsaiType1][tsaiType2] = new SigmaParameters(cEnd[tsaiType1][tsaiType2],
                        cP1[tsaiType1][tsaiType2],
                        cP2[tsaiType1][tsaiType2],
                        cValAtp1[tsaiType1][tsaiType2],
                        cValAtp2[tsaiType1][tsaiType2]);
                hbParameters[tsaiType1][tsaiType2] = new SigmaParameters(hbEnd[tsaiType1][tsaiType2],
                        hbP1[tsaiType1][tsaiType2],
                        hbP2[tsaiType1][tsaiType2],
                        hbValAtp1[tsaiType1][tsaiType2],
                        hbValAtp2[tsaiType1][tsaiType2]);
            }
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
            for (int ind1 = 0; ind1 < N_TSAI; ind1++) {
                line = br.readLine();
                if (line == null)
                    throw new RuntimeException("In " + filename + " there should be a " + N_TSAI + "x" +
                            N_TSAI + " array of doubles");
                stok = new StringTokenizer(line);
                for (int ind2 = 0; ind2 < N_TSAI; ind2++)
                    ar[ind1][ind2] = Double.valueOf(stok.nextToken().trim()).doubleValue();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }
}
