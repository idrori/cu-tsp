/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry.rotamers;

import java.util.*;
import java.io.*;

import meshi.energy.*;
import meshi.util.*;
import meshi.geometry.*;
import meshi.molecularElements.*;
import meshi.parameters.*;

/**
 * Reading and working with Dunbrack's  backbone dependent rotamer library.
 * The rotamers in this class are sorted in a decreasing order of probability.
 * <p/>
 * Note:
 * -----
 * The library file name should be given in line 19.
 * The library phi/psi resolution is given in line 20. (Dunbrack uses 10 deg)
 * See the constructor comments to learn how to create the class.
 * All the angles are in radians. Torsion angles are in the range [-PI,PI].
 * Residue types are numbered 0-19 (ALA,CYS,ASP,GLU,PHE,GLY...)
 */

public class DunbrackLib implements KeyWords {

    private String parametersFileName;
    public final double res = Math.PI / 18.0; // 10 degree bins
    public final int nbins = (int) Math.round(2 * Math.PI / res);

    private final int[] maxchi = {-1, 1, 2, 3, 2, -1, 2, 2, 4, 2, 3, 2, 2, 3, 4, 1, 1, 1, 2, 2};
    private double[][][][] prob = new double[20][nbins][nbins][];
    private double[][][][][] rot = new double[20][nbins][nbins][][];


    /**
     * The default constructor reads up to 15 rotamers or until they cover 98% of the cases.
     */
    public DunbrackLib(CommandList commands) {
        this(commands, 0.98, 15);
    }

    /**
     * The constructors needs 2 parameters to decide how many rotamers to load from the library.
     * This is neccessary because (for example) Lysine has 81 rotamers most of whom of negligebale
     * probability.
     * The first parameter is the coverage you need from the rotamers. if you put 0.9, the most
     * probable rotamers that cover at least 90% of the cases will be available.
     * The second parameter is the maximal number of rotamer to read (regardless of their coverage).
     */
    public DunbrackLib(CommandList commands, double minprob, int maxn) {

        // get the library file name
        Command command = commands.firstWord(PARAMETERS_DIRECTORY);
        String parametersFileName = command.secondWord();
        command = commands.firstWord(ROTAMER_LIBRARY);
        parametersFileName = parametersFileName + "/" + command.secondWord();

        int resnum, prevresnum = -999, ind1, ind2, cc = 0, tc, tc1, tc2;
        double phi, psi, prevphi = -999, prevpsi = -999, cummsum;
        double[] pr = new double[100];
        double[][] chi = new double[4][100];
        String line = "";
        StringTokenizer st;

        try {
            FileReader fr = new FileReader(parametersFileName);
            BufferedReader br = new BufferedReader(fr);
            System.out.println("Opened the backbone dependent library: " + parametersFileName);
            line = br.readLine();
            while (line != null) {
                st = new StringTokenizer(line);
                resnum = ResidueType.type(st.nextToken().trim()).ordinal();
                phi = Double.valueOf(st.nextToken().trim()).doubleValue();
                psi = Double.valueOf(st.nextToken().trim()).doubleValue();
                if ((psi != prevpsi) && (prevpsi > -900)) {
                    if ((prevphi != 180.0) && (prevpsi != 180.0)) {
                        ind1 = findInd(prevphi * Math.PI / 180.0);
                        ind2 = findInd(prevpsi * Math.PI / 180.0);
                        cummsum = 0;
                        for (tc = 0; (tc < cc) && (cummsum < minprob); tc++)
                            cummsum += pr[tc];
                        if (tc > maxn)
                            tc = maxn;
                        prob[prevresnum][ind1][ind2] = new double[tc];
                        rot[prevresnum][ind1][ind2] = new double[tc][maxchi[prevresnum]];
                        for (tc1 = 0; tc1 < tc; tc1++) {
                            prob[prevresnum][ind1][ind2][tc1] = pr[tc1];
                            for (tc2 = 0; tc2 < maxchi[prevresnum]; tc2++)
                                rot[prevresnum][ind1][ind2][tc1][tc2] = Math.PI * chi[tc2][tc1] / 180.0;
                        }
                    }
                    cc = 0;
                }
                st.nextToken();
                st.nextToken();
                st.nextToken();
                st.nextToken();
                st.nextToken();
                pr[cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[0][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[1][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[2][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[3][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                prevphi = phi;
                prevpsi = psi;
                prevresnum = resnum;
                cc++;
                line = br.readLine();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
        System.out.println("Rotamer library: OK");
    } // end of constructor

    public DunbrackLib(String parametersFileName) {

        // get the library file name

        this.parametersFileName = parametersFileName;
        double minprob = 0.98;
        int maxn = 15;
        int resnum, prevresnum = -999, ind1, ind2, cc = 0, tc, tc1, tc2;
        double phi, psi, prevphi = -999, prevpsi = -999, cummsum;
        double[] pr = new double[100];
        double[][] chi = new double[4][100];
        String line = "";
        StringTokenizer st;

        try {
            FileReader fr = new FileReader(parametersFileName);
            BufferedReader br = new BufferedReader(fr);
            System.out.println("Opened the backbone dependent library: " + parametersFileName);
            line = br.readLine();
            while (line != null) {
                st = new StringTokenizer(line);
                resnum = ResidueType.type(st.nextToken().trim()).ordinal();
                phi = Double.valueOf(st.nextToken().trim()).doubleValue();
                psi = Double.valueOf(st.nextToken().trim()).doubleValue();
                if ((psi != prevpsi) && (prevpsi > -900)) {
                    if ((prevphi != 180.0) && (prevpsi != 180.0)) {
                        ind1 = findInd(prevphi * Math.PI / 180.0);
                        ind2 = findInd(prevpsi * Math.PI / 180.0);
                        cummsum = 0;
                        for (tc = 0; (tc < cc) && (cummsum < minprob); tc++)
                            cummsum += pr[tc];
                        if (tc > maxn)
                            tc = maxn;
                        prob[prevresnum][ind1][ind2] = new double[tc];
                        rot[prevresnum][ind1][ind2] = new double[tc][maxchi[prevresnum]];
                        for (tc1 = 0; tc1 < tc; tc1++) {
                            prob[prevresnum][ind1][ind2][tc1] = pr[tc1];
                            for (tc2 = 0; tc2 < maxchi[prevresnum]; tc2++)
                                rot[prevresnum][ind1][ind2][tc1][tc2] = Math.PI * chi[tc2][tc1] / 180.0;
                        }
                    }
                    cc = 0;
                }
                st.nextToken();
                st.nextToken();
                st.nextToken();
                st.nextToken();
                st.nextToken();
                pr[cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[0][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[1][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[2][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                chi[3][cc] = Double.valueOf(st.nextToken().trim()).doubleValue();
                prevphi = phi;
                prevpsi = psi;
                prevresnum = resnum;
                cc++;
                line = br.readLine();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
        System.out.println("Rotamer library: OK");
    } // end of constructor


    private int findInd(double phi) {
        int result = (int) Math.round((phi + Math.PI) / res);
        if (result == nbins)
            result = 0;
        return result;
    }


    /**
     *---------------
     *Useful Methods:
     *---------------
     **/

    /**
     * Gives the number of rotamers available for a specific residue type (parameter #1) at  a certain
     * phi/psi (parameters #2,3)
     */
    public int getRotamerNum(int resnum, double phi, double psi) {
        return rot[resnum][findInd(phi)][findInd(psi)].length;
    }


    /**
     * Gives the number of Chi torsion for a specific residue type (parameter #1). For example LYS is 4.
     * ASP is 2.
     */
    public int getChiMax(int resnum) {
        return rot[resnum][0][0][0].length;
    }

    /**
     * Gives the coverage of the rotamers available for a specific residue type (parameter #1) at a certain
     * phi/psi (parameters #2,3). If for example the rotamers cover only 94% of the cases the result will be
     * 0.94
     */
    public double getCummSum(int resnum, double phi, double psi) {
        double cummsum = 0.0;
        for (int c = 0; c < getRotamerNum(resnum, phi, psi); c++)
            cummsum += prob[resnum][findInd(phi)][findInd(psi)][c];
        return cummsum;
    }

    /**
     * Gives the array of a rotamer Chi angle for a specific residue type (parameter #1) at a certain
     * phi/psi (parameters #2,3). The rotamer index is given as parameter #4. The rotamers in a certain
     * phi/psi bin are sorted in decreasing order of probabilty, and are indexed starting from 0 (for the
     * most probable one).
     * For example getRotamer(2,-140/180.0*Math.PI,130/180.0*Math.PI,1) will return an array of two doubles
     * of the Chi1,Chi2 of ASP at that phi/psi location for the second most probable rotamer.
     */
    public double[] getRotamer(int resnum, double phi, double psi, int rotNum) {
        return rot[resnum][findInd(phi)][findInd(psi)][rotNum];
    }

    /**
     * Gives the probablity for a certain rotamer (whose index is given as parameter #4) for a specific residue
     * type (parameter #1) at a certain phi/psi (parameters #2,3). The rotamers in a certain
     * phi/psi bin are sorted in decreasing order of probabilty, and are indexed starting from 0 (for the
     * most probable one).
     */
    public double getRotamerProb(int resnum, double phi, double psi, int rotNum) {
        return prob[resnum][findInd(phi)][findInd(psi)][rotNum];
    }

}
