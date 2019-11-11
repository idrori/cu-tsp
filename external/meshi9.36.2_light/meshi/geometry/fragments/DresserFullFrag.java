/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry.fragments;

import meshi.geometry.Coordinates;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.ResidueType;
import meshi.util.KeyWords;
import meshi.util.Rms;
import meshi.util.overlap.Overlap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;


public class DresserFullFrag implements KeyWords {

    private double[][][] N;
    private double[][][] CA;
    private double[][][] C;
    private double[][][] O;
    private double[][][] CB;
    private double[][][] H;
    private double[] endDis;
    private double[][] co, co2;
    public final int SEGMENT;
    public final int MAX_NUM_FRAGS;
    private final MolecularSystem molecularSystem;


    public DresserFullFrag(String fragmentFileName) {
        String line = "";
        StringTokenizer st;
        molecularSystem = new MolecularSystem();
        try {
            FileReader fr = new FileReader(fragmentFileName);
            BufferedReader br = new BufferedReader(fr);
            System.out.println("Opened the fragment library: " + fragmentFileName);
            line = br.readLine();
            st = new StringTokenizer(line);
            SEGMENT = Integer.valueOf(st.nextToken().trim()).intValue();
            if (SEGMENT / 2.0 == SEGMENT / 2)
                throw new RuntimeException("The segment length must NOT be even.");
            MAX_NUM_FRAGS = Integer.valueOf(st.nextToken().trim()).intValue();
            N = new double[MAX_NUM_FRAGS][SEGMENT][3];
            CA = new double[MAX_NUM_FRAGS][SEGMENT][3];
            C = new double[MAX_NUM_FRAGS][SEGMENT][3];
            O = new double[MAX_NUM_FRAGS][SEGMENT][3];
            CB = new double[MAX_NUM_FRAGS][SEGMENT][3];
            H = new double[MAX_NUM_FRAGS][SEGMENT][3];
            endDis = new double[MAX_NUM_FRAGS];
            co = new double[3][SEGMENT * 6];
            co2 = new double[3][SEGMENT * 6];
            for (int c = 0; c < MAX_NUM_FRAGS; c++) {
                line = br.readLine();
                st = new StringTokenizer(line);
                for (int d = 0; d < SEGMENT; d++) {
                    CA[c][d][0] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    CA[c][d][1] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    CA[c][d][2] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    N[c][d][0] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    N[c][d][1] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    N[c][d][2] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    H[c][d][0] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    H[c][d][1] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    H[c][d][2] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    CB[c][d][0] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    CB[c][d][1] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    CB[c][d][2] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    C[c][d][0] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    C[c][d][1] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    C[c][d][2] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    O[c][d][0] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    O[c][d][1] = Double.valueOf(st.nextToken().trim()).doubleValue();
                    O[c][d][2] = Double.valueOf(st.nextToken().trim()).doubleValue();
                }
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
        System.out.println("Fragment library: OK");
        for (int c = 0; c < MAX_NUM_FRAGS; c++)
            endDis[c] = (CA[c][0][0] - CA[c][SEGMENT - 1][0]) * (CA[c][0][0] - CA[c][SEGMENT - 1][0]) +
                    (CA[c][0][1] - CA[c][SEGMENT - 1][1]) * (CA[c][0][1] - CA[c][SEGMENT - 1][1]) +
                    (CA[c][0][2] - CA[c][SEGMENT - 1][2]) * (CA[c][0][2] - CA[c][SEGMENT - 1][2]);
    } // of constructor


    public AtomList dressFrag(AtomList al) {
        AtomList out = new AtomList(molecularSystem);
        int[] aux = new int[SEGMENT];
        double dis, tmpdis;
        int c;
        Atom atom;
        Rms rms;
        ResidueType aaType;
        int residueNumber;

        if (al.size() != SEGMENT)
            throw new RuntimeException("Segment length does not match.");
        for (c = 0; c < SEGMENT; c++) {
            if (!al.atomAt(c).name().equals("CA"))
                throw new RuntimeException("Only CA getAtoms can be dressed.");
            if ((c > 0) && (al.atomAt(c).residueNumber() != al.atomAt(c - 1).residueNumber() + 1))
                throw new RuntimeException("CA's in the list must be consecutive residue numbers.");
        }
/*!	dis = (al.atomAt(0).x()-al.atomAt(SEGMENT-1).x())*(al.atomAt(0).x()-al.atomAt(SEGMENT-1).x()) +
	      (al.atomAt(0).y()-al.atomAt(SEGMENT-1).y())*(al.atomAt(0).y()-al.atomAt(SEGMENT-1).y()) +
	      (al.atomAt(0).z()-al.atomAt(SEGMENT-1).z())*(al.atomAt(0).z()-al.atomAt(SEGMENT-1).z());
    for (c=0 ; c<MAX_NUM_FRAGS ; c++) 
       if (dis>endDis[c])
          break;
    System.out.print(c + " ");      */

        dis = 99999999;
        int best = -1;
        for (c = 0; c < SEGMENT; c++)
            aux[c] = c;
        for (c = 0; c < SEGMENT; c++) {
            atom = al.atomAt(c);
            co[0][c] = atom.x();
            co[1][c] = atom.y();
            co[2][c] = atom.z();
        }
        for (int fc = 0; fc < MAX_NUM_FRAGS; fc++) {
            prepareRMS(fc);
            Overlap.rmsPartial(co, co2, aux);
            tmpdis = 0.0;
            for (c = 0; c < SEGMENT; c++)
                tmpdis += (co[0][c] - co2[0][c]) * (co[0][c] - co2[0][c]) +
                        (co[1][c] - co2[1][c]) * (co[1][c] - co2[1][c]) +
                        (co[2][c] - co2[2][c]) * (co[2][c] - co2[2][c]);
            tmpdis = Math.sqrt(tmpdis / SEGMENT);

            if (tmpdis < dis) {
                dis = tmpdis;
                best = fc;
            }
        }
//!	System.out.println(best + " " + dis);	

        // Now creating the new atom list
        prepareRMS(best);
        Overlap.rmsPartial(co, co2, aux);
        for (c = 0; c < SEGMENT; c++) {
            aaType = ResidueType.type(al.atomAt(c).type());
            residueNumber = al.atomAt(c).residue().number();
            ResidueIdentifier id = new ResidueIdentifier(residueNumber);
            Residue residue = new Residue(id, aaType.nameThreeLetters());
            // Adding CA
            //                                    out.add(new Atom(co[0][c],co[1][c],co[2][c],
            //                                    "CA",aaType,al.atomAt(c).residueNumber(),-1));
            out.add(new Atom("CA", residue, aaType.caType(), new Coordinates(co[0][c], co[1][c], co[2][c]), 1, new Double(0),molecularSystem));
            // Adding C
            //                                    out.add(new Atom(co2[0][SEGMENT+c],co2[1][SEGMENT+c],co2[2][SEGMENT+c],
            //                                    "C",aaType,al.atomAt(c).residueNumber(),-1));
            out.add(new Atom("C", residue, aaType.cType(), new Coordinates(co2[0][SEGMENT + c], co2[1][SEGMENT + c], co2[2][SEGMENT + c]), 1, new Double(0),molecularSystem));
            // Adding N
            //                                    out.add(new Atom(co2[0][2*SEGMENT+c],co2[1][2*SEGMENT+c],co2[2][2*SEGMENT+c],
            //                        	         "N",aaType,al.atomAt(c).residueNumber(),-1));
            out.add(new Atom("N", residue, aaType.nType(), new Coordinates(co2[0][2 * SEGMENT + c], co2[1][2 * SEGMENT + c], co2[2][2 * SEGMENT + c]), 1, new Double(0),molecularSystem));
            // Adding CB
            if (!aaType.equals(ResidueType.GLY)) {
                //                                out.add(new Atom(co2[0][3*SEGMENT+c],co2[1][3*SEGMENT+c],co2[2][3*SEGMENT+c],
                //                                "CB",aaType,al.atomAt(c).residueNumber(),-1));
                out.add(new Atom("CB", residue, aaType.cbType(), new Coordinates(co2[0][3 * SEGMENT + c], co2[1][3 * SEGMENT + c], co2[2][3 * SEGMENT + c]),1, new Double(0),molecularSystem));
            }
            // Adding O
            //                                   out.add(new Atom(co2[0][4*SEGMENT+c],co2[1][4*SEGMENT+c],co2[2][4*SEGMENT+c],
            //                                   "O",aaType,al.atomAt(c).residueNumber(),-1));
            out.add(new Atom("O", residue, aaType.oType(), new Coordinates(co2[0][4 * SEGMENT + c], co2[1][4 * SEGMENT + c], co2[2][4 * SEGMENT + c]), 1, new Double(0),molecularSystem));
            // Adding H
            if (!aaType.equals("PRO")) {
                //                               out.add(new Atom(co2[0][5*SEGMENT+c],co2[1][5*SEGMENT+c],co2[2][5*SEGMENT+c],
                //                               "H",aaType,al.atomAt(c).residueNumber(),-1));
                out.add(new Atom("H", residue, aaType.hType(), new Coordinates(co2[0][5 * SEGMENT + c], co2[1][5 * SEGMENT + c], co2[2][5 * SEGMENT + c]), 1, new Double(0),molecularSystem));
            }
        }
        return out;
    }

    public AtomList dressProt(AtomList al) {
        AtomList out = new AtomList(molecularSystem);
        AtomList caFrag;
        AtomList allFrag;
        int first = firstRes(al);
        int last = lastRes(al);
        int dd;
        boolean goodFrag;
        boolean nextGoodFrag;
        boolean prevGoodFrag = false;
        boolean notEndFrag = true;


        for (int c = 0; c < al.size(); c++)
            if (!al.atomAt(c).name().equals("CA"))
                throw new RuntimeException("Only CA  can be dressed.");
        for (int resc = first; resc < (last - SEGMENT + 2);) {
            goodFrag = true;
            for (dd = 0; (dd < SEGMENT) && goodFrag; dd++)
                if (find(al, resc + dd) == -1)
                    goodFrag = false;
            // Checking if there is a good next segment
            nextGoodFrag = true;
            for (dd = 0; (dd < SEGMENT) && nextGoodFrag; dd++)
                if (find(al, resc + dd + SEGMENT - 2) == -1)
                    nextGoodFrag = false;
            // Dressing the fragment
            if (goodFrag) {
                caFrag = new AtomList(molecularSystem);
                for (int c = 0; c < SEGMENT; c++)
                    caFrag.add(al.atomAt(find(al, resc + c)));
                allFrag = dressFrag(caFrag);
                for (int c = 0; c < allFrag.size(); c++) {
                    if (!prevGoodFrag && (allFrag.atomAt(c).residueNumber() == resc))
                        out.add(allFrag.atomAt(c));
                    if (notEndFrag && (allFrag.atomAt(c).residueNumber() > resc) && (allFrag.atomAt(c).residueNumber() < resc + SEGMENT - 1))
                        out.add(allFrag.atomAt(c));
                    if (!notEndFrag && (find(out, allFrag.atomAt(c).name, allFrag.atomAt(c).residueNumber()) == -1))
                        out.add(allFrag.atomAt(c));
                }
                prevGoodFrag = true;
            } else
                prevGoodFrag = false;

            if (goodFrag && nextGoodFrag) {
                resc += SEGMENT - 2;
                notEndFrag = true;
            }
            if (goodFrag && !nextGoodFrag && notEndFrag) {
                resc += dd - 3;
                notEndFrag = false;
            } else if (goodFrag && !nextGoodFrag && !notEndFrag) {
                resc++;
                notEndFrag = true;
            }
            if (!goodFrag) {
                resc++;
                notEndFrag = true;
            }
        }  // of checking the chain

        return out;
    }

    private void prepareRMS(int fc) {
        int c;
        for (c = 0; c < SEGMENT; c++) {
            co2[0][c] = CA[fc][c][0];
            co2[1][c] = CA[fc][c][1];
            co2[2][c] = CA[fc][c][2];
        }
        for (c = 0; c < SEGMENT; c++) {
            co2[0][SEGMENT + c] = C[fc][c][0];
            co2[1][SEGMENT + c] = C[fc][c][1];
            co2[2][SEGMENT + c] = C[fc][c][2];
        }
        for (c = 0; c < SEGMENT; c++) {
            co2[0][2 * SEGMENT + c] = N[fc][c][0];
            co2[1][2 * SEGMENT + c] = N[fc][c][1];
            co2[2][2 * SEGMENT + c] = N[fc][c][2];
        }
        for (c = 0; c < SEGMENT; c++) {
            co2[0][3 * SEGMENT + c] = CB[fc][c][0];
            co2[1][3 * SEGMENT + c] = CB[fc][c][1];
            co2[2][3 * SEGMENT + c] = CB[fc][c][2];
        }
        for (c = 0; c < SEGMENT; c++) {
            co2[0][4 * SEGMENT + c] = O[fc][c][0];
            co2[1][4 * SEGMENT + c] = O[fc][c][1];
            co2[2][4 * SEGMENT + c] = O[fc][c][2];
        }
        for (c = 0; c < SEGMENT; c++) {
            co2[0][5 * SEGMENT + c] = H[fc][c][0];
            co2[1][5 * SEGMENT + c] = H[fc][c][1];
            co2[2][5 * SEGMENT + c] = H[fc][c][2];
        }
    }

    private int firstRes(AtomList al) {
        int result = 100000000;
        for (int c = 0; c < al.size(); c++)
            if (al.atomAt(c).residueNumber() < result)
                result = al.atomAt(c).residueNumber();
        return result;
    }

    private int lastRes(AtomList al) {
        int result = -100000000;
        for (int c = 0; c < al.size(); c++)
            if (al.atomAt(c).residueNumber() > result)
                result = al.atomAt(c).residueNumber();
        return result;
    }

    private int find(AtomList al, int res) {
        for (int c = 0; c < al.size(); c++)
            if (al.atomAt(c).residueNumber() == res)
                return c;
        return -1;
    }

    private int find(AtomList al, String name, int res) {
        for (int c = 0; c < al.size(); c++)
            if (al.atomAt(c).name().equals(name) && (al.atomAt(c).residueNumber() == res))
                return c;
        return -1;
    }

} // of class

