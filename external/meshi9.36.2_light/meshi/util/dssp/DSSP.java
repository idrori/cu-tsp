/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.dssp;

import meshi.parameters.ResidueType;
import meshi.sequences.*;
import meshi.util.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

public class DSSP {
    protected String fileName;
    protected char[] aa;
    protected char[] chain;
    protected char[] chains;
    protected int[] resNum;
    protected char[] ss;
    protected  double[] solvACC;
    protected  double[] fullSol = {115, 135, 150, 190, 210,
            75, 195, 175, 200, 170,
            185, 160, 145, 180, 225,
            115, 140, 155, 255, 230};
    protected double hbPrecentage;
    protected double parallelPrecentage;
    protected double antyparallelPrecentage;

    public void printPresentage() {
        System.out.println("Total hb: " + hbPrecentage);
        System.out.println("parallel: " + parallelPrecentage);
        System.out.println("AntiParallel: " + antyparallelPrecentage);

        System.out.println("Number Of Betta Residues: ");
    }

    public DSSP(String dsspFileName) {
        this.fileName = dsspFileName;
        readDSSP();
    }

    /*
    *This function returns true if the DSSP assignment of residue number 'res'
    *appears in the given char list. Else false is returned.
    */
    public boolean inList(int res, char[] SSlist) {
        char ch;
        int i;

        ch = SSofRes(res);
        for (i = 0; i < SSlist.length; i++) {
            if (ch == SSlist[i])
                return true;
        }
        return false;
    }


    /*
    *This function returns true if the residue given is not in the edge of a structure.
    *i.e. the residues +/-1 has the same SS as it has.
    */
    public boolean notInEdge(int res) {
        char chp1 = ' ', chm1 = ' ', ch = ' ';
        int i;
        for (i = 0; (i < resNum.length); i++) {
            if (resNum[i] == (res - 1))
                chm1 = ss[i];
            if (resNum[i] == (res))
                ch = ss[i];
            if (resNum[i] == (res + 1))
                chp1 = ss[i];
        }
        return (chp1 == ch) && (ch == chm1);
    }

    /*
    *Returns the SS structure of residue number 'res'.
    *If this residue number doesn't exist - 'X' is returned.
    */
    public char SSofRes(int res) {
        int i;
        for (i = 0; i < resNum.length; i++) {
            if (resNum[i] == res)
                return ss[i];
        }
        return 'X';
    }


    /*
    *Returns the AA of residue number 'res'.
    *If this residue number doesn't exist - 'X' is returned.
    */
    public char AAofRes(int res) {
        int i;
        for (i = 0; i < resNum.length; i++) {
            if (resNum[i] == res)
                return aa[i];
        }
        return 'X';
    }

    /*
    *Returns the solvation accessibilty of residue number 'res'.
    *If this residue number doesn't exist - (-1) is returned.
    */
    public double ACCofRes(int res) {
        int i;
        for (i = 0; i < resNum.length; i++) {
            if (resNum[i] == res)
                return solvACC[i];
        }
        return -1;
    }


    /*
    *Returns the relative solvation accessibilty of residue number 'res'.
    *If this residue number doesn't exist - (-1) is returned.
    */
    public double relACCofRes(int res) {
        int i;
        for (i = 0; i < resNum.length; i++) {
            if (resNum[i] == res)
                return solvACC[i] / fullSol[ResidueType.type(aa[i]).ordinal()];
        }
        return -1;
    }

    public String getACC() {
        String ans = "";
        for (int i = 0; i < resNum.length; i++)
            if (resNum[i] > -999)
                ans = ans + solvACC[i] + ";";
        return ans;
    }

    public String[] getRelativeACC(char chainChar) {
        ArrayList<String> out = new ArrayList();
        for (int iResidue = 0; iResidue < resNum.length; iResidue++)
            if (chain[iResidue] == chainChar) {
                if (resNum[iResidue] > -999) {
                    if ((aa[iResidue] == 'a') || (aa[iResidue] == 'b') || (aa[iResidue] == 'c') || (aa[iResidue] == 'd'))
                        aa[iResidue] = 'C';// Cysteins involved in disulphide bonds
                    if ((ResidueType.type(aa[iResidue]).ordinal() > 19) || (ResidueType.type(aa[iResidue]).ordinal() < 0))
                        throw new RuntimeException("weird residue " + aa[iResidue] + " type number " + ResidueType.type(aa[iResidue]).ordinal());
                    if ((int) ((solvACC[iResidue] / fullSol[ResidueType.type(aa[iResidue]).ordinal()]) + 0.5) >= 1)
                        out.add("A" + (solvACC[iResidue] / ResidueType.type(aa[iResidue]).maxExposedSurfaceArea));
                    else
                        out.add("B" + (solvACC[iResidue] / ResidueType.type(aa[iResidue]).maxExposedSurfaceArea));
                    if (fullSol[ResidueType.type(aa[iResidue]).ordinal()] != ResidueType.type(aa[iResidue]).maxExposedSurfaceArea)
                        throw new RuntimeException("This is weird ");
                }
            }
        return out.toArray(new String[1]);
    }

    public String getAA(char chainChar) {
        String ans = "";
        for (int iResidue = 0; iResidue < resNum.length; iResidue++) {
            if (chain[iResidue] == chainChar) {
                if ((aa[iResidue] >= 'a') && (aa[iResidue] <= 'z'))
                    aa[iResidue] = 'C';// Cysteins involved in disulphide bonds
                if (resNum[iResidue] > -999)
                    ans = ans + aa[iResidue];
                else
                    throw new RuntimeException("res " + iResidue + " is -999" + "\n" +
                            ans);
            }
        }
        return ans;
    }

    public void printSS() {
        for (int i = 0; i < resNum.length; i++) {
            if (resNum[i] > -999)
                System.out.println(ss[i]);
            //System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
        }
    }

    public String getSSInOneLine() {
        String ans = "";
        for (int i = 0; i < resNum.length; i++) {
            if (resNum[i] > -999)
                ans = ans + ss[i];
            /*if(ss[i] == 'B' | ss[i] =='E')
             ans = ans+'E';
             else if(ss[i] == 'H')
             ans = ans+'H';
             else if(ss[i] == 'S' | ss[i] == 'G' | ss[i] == 'I' | ss[i] == 'T')
             ans = ans+'A';
             else
             ans = ans+'C';*/

            //System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
        }
        return ans;
    }

    public String getSSInOneLineConventional(char chainChar) {
        String ans = "";
        for (int iResidue = 0; iResidue < resNum.length; iResidue++) {
            if (chain[iResidue] == chainChar) {
                if (resNum[iResidue] > -999)

                    if (ss[iResidue] == 'B' | ss[iResidue] == 'E')
                        ans = ans + 'E';
                    else if (ss[iResidue] == 'H'| ss[iResidue] == 'G' | ss[iResidue] == 'I' )
                        ans = ans + 'H';
                        //else if(ss[iResidue] == 'S' | ss[i] == 'T')
                        //  ans = ans+'A';//TODO - add this !
                    else
                        ans = ans + 'C';

                //System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
            }
        }
        return ans;
    }

    public int getNumberOfBetaResidues() {
        int ans = 0;
        for (int i = 0; i < resNum.length; i++) {
            if (resNum[i] > -999 && ss[i] == 'E')  //TODO(ss[i] == 'B' |
                ans++;
        }
        return ans;
    }

    protected static   String getDsspFileName(String fileName) {
        String testedFileNames = "";
        File dsspFile = null;
        while ((dsspFile == null) || (!dsspFile.exists())) {
             testedFileNames += " "+fileName+" ";
             dsspFile = new File(fileName+".dssp");
            if (!dsspFile.exists()) {
                if (fileName.lastIndexOf('.')<0) {
                    throw new RuntimeException("Failed to find dssp file in any of "+testedFileNames);
                }
                fileName = fileName.substring(0,fileName.lastIndexOf('.'));
            }
        }
         return dsspFile.getAbsolutePath();
    }

    public int length() {
        return aa.length;
    }
    public void readDSSP() {
        int nResidues = 0;
        int nChains = 0;
        char currentChain = '$';
        boolean cont = true;
        String line;
        try {

            String dsspFileName = getDsspFileName(fileName);
            FileReader fr = new FileReader(dsspFileName);
            BufferedReader fdssp = new BufferedReader(fr);
            line = fdssp.readLine();
            while ((line != null) && (cont)) {
                if (line.length() > 25)
                    if (line.substring(0, 25).compareTo("  #  RESIDUE AA STRUCTURE") == 0) {
                        cont = false;
                    } else {
                        String word = line.substring(50, 65);
                        if (word.equalsIgnoreCase(" Parallel BRIDG"))
                            parallelPrecentage = Double.parseDouble(line.substring(6, 10));
                        else if (word.equalsIgnoreCase("iParallel BRIDG"))
                            antyparallelPrecentage = Double.parseDouble(line.substring(6, 10));
                        else if (word.equalsIgnoreCase("E O(I)-->H-N(J)"))
                            hbPrecentage = Double.parseDouble(line.substring(6, 10));
                    }
                line = fdssp.readLine();
            }
            while (line != null) {
                if ((line != null) && (line.indexOf("!")) > -1) {
                    line = fdssp.readLine();
                    continue;
                }

                if (line.charAt(11) != currentChain) {
                    nChains++;
                    currentChain = line.charAt(11);
                }
                nResidues++;
                line = fdssp.readLine();
            }
            fdssp.close();
            currentChain = '$';
            fr = new FileReader(dsspFileName);
            fdssp = new BufferedReader(fr);
            resNum = new int[nResidues];
            ss = new char[nResidues];
            aa = new char[nResidues];
            chains =  new char[nChains];
            chain = new char[nResidues];
            solvACC = new double[nResidues];
            for (int cc = 0; cc < resNum.length; cc++)
                resNum[cc] = -999;
            int iResidue = 0;
            int iChain = 0;
            cont = true;
            line = fdssp.readLine();
            while ((line != null) && (cont)) {
                if (line.length() > 25)
                    if (line.substring(0, 25).compareTo("  #  RESIDUE AA STRUCTURE") == 0) {
                        cont = false;
                    }
                line = fdssp.readLine();
            }
            while (line != null) {
                try {
                    resNum[iResidue] = (Integer.valueOf(line.substring(5, 10).trim()));
                    chain[iResidue] = line.charAt(11);
                    if (chain[iResidue] != currentChain) {
                        currentChain   = chain[iResidue];
                        chains[iChain] = currentChain;
                        iChain++;
                    }
                    aa[iResidue] = line.charAt(13);
                    ss[iResidue] = line.charAt(16);
                    if (ss[iResidue] == ' ')
                        ss[iResidue] = 'C';
                    solvACC[iResidue] = (Double.valueOf(line.substring(35, 38).trim()));

                }
                catch (Exception e) {
                    iResidue--;
                }
                iResidue++;
                line = fdssp.readLine();
            }
            fdssp.close();
        } // of the try
        catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        } // of the catch
    }

    public SequenceList[] sequenceLists() {
        SequenceList[] out = new SequenceList[chains.length];
        for (int iChain = 0; iChain < chains.length; iChain++) {
            SequenceList seqList = new SequenceList();
            MeshiSequence aaSeq = new ResidueSequence(getAA(chains[iChain]), "aa sequence of " + fileName);
            seqList.add(aaSeq);
            MeshiSequence ssSeq = new SecondaryStructureSequence(getSSInOneLineConventional(chains[iChain]), "ss sequence of " + fileName);
            seqList.add(ssSeq);
            String[] tempAccSeq = getRelativeACC(chains[iChain]);
            MeshiSequence accSeq = new AccesibilitySequence(tempAccSeq, "accesibility sequence of " + fileName);
            seqList.add(accSeq);
            out[iChain] = seqList;
        }
        return out;
    }
}


/*
 * The dssp line:
 1         2         3         4         5         6         7         8         9         0
 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
 18   23 A I  T <4 S+     0   0   64     -3,-3.4    -2,-0.2    -4,-0.2     7,-0.2   0.722 140.9  39.8-114.4 -49.1   11.0
   
*/
