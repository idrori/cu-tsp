/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

/**
 * Converts Residue names from three letter to one letter format.
 */
public class ResidueNameConvertor {
    public static char three2one(String threeLName) {
        char ans = '*';
        threeLName = threeLName.trim();
        if (threeLName.length() == 3) {
            if (threeLName.equals("ALA"))
                ans = 'A';
            else {
                if (threeLName.equals("CYS"))
                    ans = 'C';
                else {
                    if (threeLName.equals("ASP"))
                        ans = 'D';
                    else {
                        if (threeLName.equals("GLU"))
                            ans = 'E';
                        else {
                            if (threeLName.equals("PHE"))
                                ans = 'F';
                            else {
                                if (threeLName.equals("GLY"))
                                    ans = 'G';
                                else {
                                    if (threeLName.equals("HIS"))
                                        ans = 'H';
                                    else {
                                        if (threeLName.equals("ILE"))
                                            ans = 'I';
                                        else {
                                            if (threeLName.equals("LYS"))
                                                ans = 'K';
                                            else {
                                                if (threeLName.equals("LYS"))
                                                    ans = 'K';
                                                else {
                                                    if (threeLName.equals("LEU"))
                                                        ans = 'L';
                                                    else {
                                                        if (threeLName.equals("MET"))
                                                            ans = 'M';
                                                        else {
                                                            if (threeLName.equals("ASN"))
                                                                ans = 'N';
                                                            else {
                                                                if (threeLName.equals("PRO"))
                                                                    ans = 'P';
                                                                else {
                                                                    if (threeLName.equals("GLN"))
                                                                        ans = 'Q';
                                                                    else {
                                                                        if (threeLName.equals("ARG"))
                                                                            ans = 'R';
                                                                        else {
                                                                            if (threeLName.equals("SER"))
                                                                                ans = 'S';
                                                                            else {
                                                                                if (threeLName.equals("THR"))
                                                                                    ans = 'T';
                                                                                else {
                                                                                    if (threeLName.equals("VAL"))
                                                                                        ans = 'V';
                                                                                    else {
                                                                                        if (threeLName.equals("TRP"))
                                                                                            ans = 'W';
                                                                                        else {
                                                                                            if (threeLName.equals("TYR"))
                                                                                                ans = 'Y';
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return ans;
    }
}
