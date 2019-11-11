package meshi.util.dssp;

import meshi.parameters.LocalStructure;
import meshi.sequences.*;

import java.io.BufferedReader;
import java.io.FileReader;

/**
 * Created by user on 14/02/2017.
 */
public class DSSPFull extends DSSP {

    protected LocalStructure[] dls;

    public DSSPFull(String dsspFileName) {
        super(dsspFileName);
    }

    public String getFullSSOfChain(char chainChar) {
        String ans = "";
        for (int iResidue = 0; iResidue < resNum.length; iResidue++) {
            if (chain[iResidue] == chainChar) {
                if (resNum[iResidue] > -999)
                    ans += String.format("%4s %8s \n",resNum[iResidue],dls[iResidue].toString());
            }
        }
        return ans;
    }


    public MeshiSequence[] getResidueSequenceWithFullSs() {
        MeshiSequence[] out = new MeshiSequence[chains.length];
        for (int iChain = 0; iChain < chains.length; iChain++) {
            MeshiSequence aaSeq = new LocalStructureSequence(getAA(chains[iChain]),getFullSs(chains[iChain]), "aa sequence of " + fileName);
            out[iChain] = aaSeq;
        }
        return out;
    }
    @Override
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
            MeshiSequence ssFullSeq = new LocalStructureSequence(getAA(chains[iChain]),getFullSs(chains[iChain]), "full ss sequence of " + fileName);
            seqList.add(ssFullSeq);
            out[iChain] = seqList;
        }
        return out;
    }

    public LocalStructure[] getFullSs(char chainChar) {
        String chainSequence = getAA(chainChar);
        LocalStructure[] ans = new LocalStructure[chainSequence.length()];
        int iSs = 0;
        for (int iResidue = 0; iResidue < resNum.length; iResidue++) {
            if (chain[iResidue] == chainChar && resNum[iResidue] > -999) {
                ans[iSs] = dls[iResidue];
                iSs++;
            }
        }
        return ans;
    }


    @Override
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
            };
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
            dls = new LocalStructure[nResidues];
            char[] dlsRawStructure;
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
                    dlsRawStructure = getResidueFullDsspStructure(line);
                    dls[iResidue] = new LocalStructure(dlsRawStructure);
                    if (ss[iResidue] == ' ')
                        ss[iResidue] = 'C';
                    solvACC[iResidue] = (Double.valueOf(line.substring(35, 38).trim()));

                }
                catch (Exception e) {
                    System.out.println(e.getMessage());
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
    private char[] getResidueFullDsspStructure(String dsspResidueLine){
        String str = dsspResidueLine.substring(16,17) + dsspResidueLine.substring(18,25);
        str = str.replace(' ','_');
        return str.toCharArray();
    }

}
