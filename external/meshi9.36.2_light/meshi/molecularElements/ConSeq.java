package meshi.molecularElements;

import meshi.util.file.MeshiLineReader;

import java.io.File;
import java.io.IOException;

/**
 */
public class ConSeq extends Conservation {
    public ConSeq(double conservationScore) {
        super(conservationScore);
    }
    public double conservationWeight() {
    //   if (2.5-conservationScore< 0) throw new RuntimeException("This is weird "+conservationScore);
        return Math.exp(-10*conservationScore);
    }

    public static boolean setConSeq(Chain chain, String pdbFileName ) throws IOException{

        String fileName =   getConSeqFile(pdbFileName);
        if (fileName == null) {
            System.out.println("failed to find conSeq file.");
            return false;
        }
        MeshiLineReader conSeqReader = new MeshiLineReader(fileName);

        String line = "";
        while (line.indexOf("(normalized)")<0) {
            line = conSeqReader.readLine();
        }
        boolean done = false;
        double maxScore = -1000;
        while(!done) {
            line = conSeqReader.readLine().replace('\t',' ');
              String[] words = line.split(" ");
            int i = 0;
            String resNameNumber = "";
            String scoreString = "";
            for( String s :words) {
                if (s.length() > 0) {
                    i++;
                    if (i == 4) scoreString = s;
                    if (i == 3) resNameNumber = s;
                }
            }
            if ((scoreString.length() == 0)|| (resNameNumber.length() == 0))
                done = true;
            else {
                String resName = resNameNumber.substring(0,3);
                String resNumberAndChain = resNameNumber.substring(3,resNameNumber.length());
                String[]  resNumberAndChainSplitted = resNumberAndChain.split(":");
                int resNumber = Integer.valueOf(resNumberAndChainSplitted[0]);
                double score = Double.valueOf(scoreString);
		if (score > maxScore) maxScore = score;
		if (resNumber < chain.size()) {
                	Residue residue = chain.residueAt(resNumber);
			if (! residue.dummy()){
                		if (!residue.name.equals(resName)) throw new RuntimeException("line: "+line+"\nResidue: "+residue+"\nresNumber: "+resNumber);
                		residue.setConservation(new ConSeq(score));
			}
		}
            }
        }
	for (Residue residue : chain) {
		if (residue.conservation() == null)
			residue.setConservation(new ConSeq(maxScore));
	}
        return true;
    }

     private static   String getConSeqFile(String fileName) {
        String testedFileNames = "";
        File conSeqFile = null;
        while ((conSeqFile == null) || (!conSeqFile.exists())) {
             testedFileNames += " "+fileName+" ";
             conSeqFile = new File(fileName+".conSeq");
            if (!conSeqFile.exists()) {
                if (fileName.lastIndexOf('.')<0) {
                    System.out.println("Failed to find conSeq file in any of "+testedFileNames);
                    return null;
                }
                fileName = fileName.substring(0,fileName.lastIndexOf('.'));
            }
        }
         return conSeqFile.getAbsolutePath();
    }
}
