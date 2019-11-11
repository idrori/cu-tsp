/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB;

import meshi.util.*;
import meshi.parameters.*;

import java.util.Formatter;

/**
 * from http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
 * <p/>
 * COLUMNS        DATA TYPE       FIELD         DEFINITION
 * ---------------------------------------------------------------------------------
 * 1 -  6        Record name     "ATOM  "
 * <p/>
 * 7 - 11        Integer         serial        Atom serial number.
 * <p/>
 * 13 - 16        Atom            name          Atom name.
 * <p/>
 * 17             Character       altLoc        Alternate location indicator.
 * <p/>
 * 18 - 20        Residue name    resName       Residue name.
 * <p/>
 * 22             Character       chainID       Chain identifier.
 * <p/>
 * 23 - 26        Integer         resSeq        Residue sequence number.
 * <p/>
 * 27             AChar           iCode         Code for insertion of residues.
 * <p/>
 * 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
 * Angstroms.
 * <p/>
 * 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
 * Angstroms.
 * <p/>
 * 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
 * Angstroms.
 * <p/>
 * 55 - 60        Real(6.2)       occupancy     Occupancy.
 * <p/>
 * 61 - 66        Real(6.2)       tempFactor    Temperature factor.
 * <p/>
 * 73 - 76        LString(4)      segID         Segment identifier, left-justified.
 * <p/>
 * 77 - 78        LString(2)      element       Element symbol, right-justified.
 * <p/>
 * 79 - 80        LString(2)      charge        Charge isOn the atom.
 */
public class PdbLine {
    private String line;
    private int number;
    private String name;
    private String alternateLocation;
    private String residueName;
    private String chain;
    private int residueNumber;
    private double x;
    private double y;
    private double z;
    private double occupancy;
    private double temperatureFactor;


    public PdbLine(String line) {
        this.line = line;
        if (isAnAtomOrHeteroAtom()) {
            try {
                number = Integer.valueOf(line.substring(6, 11).trim()).intValue();
                name = line.substring(12, 16).trim();
                alternateLocation = line.substring(16, 17).trim();
                residueName = line.substring(17, 20).trim();
                residueNumber = Integer.valueOf(line.substring(22, 26).trim());
                chain = line.substring(21, 22);
                x = Double.valueOf(line.substring(30, 38).trim()).doubleValue();
                y = Double.valueOf(line.substring(38, 46).trim()).doubleValue();
                z = Double.valueOf(line.substring(46, 54).trim()).doubleValue();
                if (line.length() > 65)
                    try {
                        temperatureFactor = Double.valueOf(line.substring(60, 66).trim()).doubleValue();
                    }
                    catch (NumberFormatException nfe) {
                        temperatureFactor = 0.0;
                    }
                else
                    temperatureFactor = 0.0;

                if (line.length() > 54)
                    try {
                        occupancy = Double.valueOf(line.substring(55, 61).trim()).doubleValue();
                    }
                    catch (NumberFormatException nfe) {
                        occupancy = 0.01;
                    }
                else
                    occupancy = 0.99;
            }
            catch (Exception ex) {
                throw new RuntimeException("PdbLine Error - Badly formatted PdbLine:\n" + line + "\n" + ex);
            }
        }
    }
    public boolean isTer() {
        return line.startsWith("TER");
    }
    public boolean isEnd() {
        return line.startsWith("END");
    }

    public PdbLine(int number, String name, String alternativeLocation, String residueName,
                   String chain, int residueNumber,
                   double x, double y, double z, double occupancy,
                   double temperatureFactor) {
        this.number = number;
        this.name = name;
        this.alternateLocation = alternateLocation;
        this.residueName = residueName;
        this.chain = chain;
        this.residueNumber = residueNumber;
        this.x = x;
        this.y = y;
        this.z = z;
        this.occupancy = occupancy;
        this.temperatureFactor = temperatureFactor;
        this.line = printLine(number, name, alternativeLocation, residueName,
                chain, residueNumber,
                x, y, z, occupancy,
                temperatureFactor);
    }

    public PdbLine() {
        this.number = -1;
        this.name = "XX";
        this.alternateLocation = "";
        this.residueName = "XXX";
        this.chain = "X";
        this.residueNumber = -1;
        this.x = -999.999;
        this.y = -999.999;
        this.z = -999.999;
        this.occupancy = -1;
        this.temperatureFactor = -1;
        this.line = printLine(number, name, residueName, alternateLocation,
                chain, residueNumber,
                x, y, z, occupancy,
                temperatureFactor);
    }

    // A potential bug - the method gets  alternativeLocation but does nothing with it.   

    private String printLine(int number, String name, String alternativeLocation, String residueName,
                             String chain, int residueNumber,
                             double x, double y, double z, double occupancy,
                             double temperatureFactor) {
        StringBuffer res = new StringBuffer();
        Formatter f = new Formatter(res);

        // occupancy is between 0.00 and 1.00
        if (occupancy < 0)
            occupancy = 0;
        else if (occupancy > 1)
            occupancy = 1;
        // temperature factor's threshold are -99 and 999
        if (temperatureFactor <= -99)
            temperatureFactor = 99.0;
        else if (temperatureFactor >= 999.0)
            temperatureFactor = 999.0;

        switch (name.length()) {
            case 1:
                name = " " + name + "  ";
                break;
            case 2:
                name = " " + name + " ";
                break;
            case 3:
                name = " " + name;
                break;
            case 4:
                name = name;
                break;
        }
        try {
            f.format("%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s",
                    "ATOM",
                    number,
                    name,
                    " ",
                    residueName,
                    chain,
                    residueNumber,
                    " ",
                    x, y, z,
                    occupancy, temperatureFactor,
                    " ", " ", " ");
        }
        catch (OutOfMemoryError er) {
            System.out.println("**************************\n" + "Failed to format\n" + "ATOM" + "\n" +
                    number + "\n" +
                    name + "\n" +
                    " " + "\n" +
                    residueName + "\n" +
                    chain + "\n" +
                    residueNumber + "\n" +
                    " " + "\n" +
                    x + "\n" + y + "\n" + z + "\n" +
                    occupancy + "\n" + temperatureFactor + "\n");
            throw er;
        }
        f.flush();
        f.close();
        return (line = res.toString());
    }

    public boolean isAnAtom() {
        return ((line.length() >= 54) && line.startsWith("ATOM"));
    }

    public boolean isAHeteroAtom() {
        return ((line.length() >= 54) && line.startsWith("HETATM"));
    }

    public boolean isAnAtomOrHeteroAtom() {
        return (isAnAtom() || isAHeteroAtom());
    }

    public void needsToBeAnAtom() {
        if (!isAnAtomOrHeteroAtom())
            throw new MeshiException("PdbLine error:\n" +
                    "needs to be an atom:\n" +
                    this + "\n");
    }

    public boolean isAComment() {
        return (!isAnAtomOrHeteroAtom());
    }

    public boolean isSEQRES() {
        return line.startsWith("SEQRES");
    }

    public String toString() {
        return line;
    }

    public double x() {
        needsToBeAnAtom();
        return x;
    }

    public double y() {
        needsToBeAnAtom();
        return y;
    }

    public double z() {
        needsToBeAnAtom();
        return z;
    }

    public String chain() {
        needsToBeAnAtom();
        return chain;
    }

    public String residueName() {
        needsToBeAnAtom();
        return residueName;
    }

    public String name() {
        return name;
    }

    public Integer residueNumber() {
        needsToBeAnAtom();
        return residueNumber;
    }

    public int number() {
        needsToBeAnAtom();
        return number;
    }
    /*
     *Check if this is a MODEL line.
     */

    public boolean isAModel() {
        return line.startsWith("MODEL");
    }
    /*
    *If this is a MODEL line, the function returns the model number, otherwise it returns -1.
    */

    public int getModel() {
        if (isAModel())
            return Integer.valueOf(line.substring(6).trim()).intValue();
        else
            return -1;
    }

    public double temperatureFactor() {
        needsToBeAnAtom();
        return temperatureFactor;
    }

    public double occupancy() {
        needsToBeAnAtom();
        return occupancy;

    }

    public String alternateLocation() {
        needsToBeAnAtom();
        return alternateLocation;
    }

    public AtomType type() {
        ResidueType residueType = ResidueType.type(residueName);
        AtomType atomType = AtomType.type(residueType.nameOneLetter(), name);
        return atomType;
    }

    public int length() {
        return line.length();
    }
}
