/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.twoTorsions;

import meshi.energy.Parameters;
import meshi.parameters.SecondaryStructure;

import java.io.BufferedReader;
import java.io.FileReader;

/**
 * The parameters for the potential are stored in the file that is the input to the constructor.
 * <p/>
 * The parameter file holds a long column of doubles written in ASCII format. The first two doubles are the 'limit'
 * and 'resolution'. Next there are 21 numbers that are in either 0 or 1 followed by the doubles for 'coef' (see
 * bellow what these variables mean). The doubles list of 'coef' is arranged by: Amino-acids, then by the different
 * polynomial pieces, then by the coefficients for each piece. See the code of 'readCoefficients()' for more clarifications.
 * <p/>
 * limit - is the limit of the ramachandran map for both phi and psi. It is usually PI, so that both phi and psi
 * are in the range [-PI,PI].
 * resolution - The 2D Ramachandran map is divided into equal squares with side of size "resolution".
 * 21 numbers - For the 20 amino acids + the omni-amino acid (an average of the other 20). If a number is 1 than this
 * amino acid is 'kosher'. If it is 0 than this amino acid can not be used in this energy (for example CHI1 in glycins),
 * and using it will lead to a run time error.
 * coefs[Number of 'kosher' amino acids][Number of squares][16] - is a 3D array that contains the coefficients for
 * the cubic splines. Its first index is the type of residue since each residue has its own map. The residue types
 * are numbered according to our convention.
 * The second index referes to a square in the map. The numbering starts with the square at [-limit,-limit] ,
 * continues along psi to [-limit,-limit+reolution], and so isOn until [-limit,+limit].The numbering than increase
 * the phi, and continue with square [-limit+resolution,-limit]. This continues until the last square being [limit,limit].
 * The third index enumerates the coeficients in the following sequence: 1 , phi , phi^2 , phi^3 , psi , psi*phi ,
 * psi*phi^2 , psi*phi^3 , psi^2 , psi^2*phi , psi^2*phi^2 , psi^2*phi^3 , psi^3 , psi^3*phi , psi^3*phi^2 , psi^3*phi^3.
 */

public class TwoTorsionsParameters implements Parameters {

    public String parametersFileName;
    public String torsion1Name;
    public String torsion2Name;
    public SecondaryStructure secondaryStructure;
    public double limit;
    public double resolution;
    public double[][][] coef;
    public int[] mapping = new int[21];
    private int numOfEntries = 0;


    public TwoTorsionsParameters(String parametersFileName) {
        this.parametersFileName = parametersFileName;
        int c1, c2, c3;
        String line;

        try {
            FileReader fr = new FileReader(parametersFileName);
            BufferedReader br = new BufferedReader(fr);

            line = br.readLine();
            torsion1Name = line.trim();
            line = br.readLine();
            torsion2Name = line.trim();
            line = br.readLine();
            String ss = line.trim();
            if (ss.equals("SHEET"))
                secondaryStructure = SecondaryStructure.SHEET;
            else if (ss.equals("HELIX"))
                secondaryStructure = SecondaryStructure.HELIX;
            else if (ss.equals("COIL"))
                secondaryStructure = SecondaryStructure.COIL;
            else if (ss.equals("ALL"))
                secondaryStructure = SecondaryStructure.ALL;
            else throw new RuntimeException("no match to secondary structure at TwoTorsionsParameters");

            line = br.readLine();
            limit = Double.valueOf(line.trim());
            line = br.readLine();
            resolution = Double.valueOf(line.trim());
            Long lo = 4 * Math.round(limit / resolution) * Math.round(limit / resolution);
            int Nsquares = lo.intValue();
            c2 = 0;
            for (c1 = 0; c1 < 21; c1++) {
                line = br.readLine().trim();
                if ((line.compareTo("0") != 0) && (line.compareTo("1") != 0))
                    throw new RuntimeException("Mapping values (begining of file)" +
                            parametersFileName + " need to be with values {0,1) only\n");
                c3 = Integer.valueOf(line.trim());
                if (c3 == 1) {
                    mapping[c1] = c2;
                    c2++;
                    numOfEntries++;
                } else if (c3 == 0) {
                    mapping[c1] = -1;
                } else {
                    throw new RuntimeException("Values in lines 3-24 of file:" +
                            parametersFileName + "\nneed to be with values {0,1) only\n");
                }
            }
            coef = new double[numOfEntries][Nsquares][16];
            for (c1 = 0; c1 < numOfEntries; c1++) {
                for (c2 = 0; c2 < Nsquares; c2++) {
                    for (c3 = 0; c3 < 16; c3++) {
                        line = br.readLine();
                        coef[c1][c2][c3] = Double.valueOf(line.trim());
                    }
                }
            }
            fr.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }


    public String toString() {
        return "TwoTorsionsParameters\n" +
                "\t torsion angle 1   = " + torsion1Name +
                "\t torsion angle 2   = " + torsion2Name +
                "\t secondary structure   = " + secondaryStructure;
    }


}
