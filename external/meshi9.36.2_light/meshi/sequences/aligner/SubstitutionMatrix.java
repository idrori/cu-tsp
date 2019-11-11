/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

public abstract class SubstitutionMatrix {
    public final static int asciiA=65;
    public final static int numAminoAcids=20, numWildChars=3, numEnglishLetters=26;
    protected int substitutionMatrix[][]=new int[numAminoAcids+numWildChars][numAminoAcids+numWildChars];
    protected char[] order=new char[numAminoAcids+numWildChars];
    protected int letterToIndex[]=new int[numEnglishLetters];
    protected double gapStartPenalty;
    protected double inGapPenalty;
    protected String[] matrixLineSplit;

    public SubstitutionMatrix(double gapStartPenalty, double inGapPenalty) {
        this.gapStartPenalty = gapStartPenalty;
        this.inGapPenalty = inGapPenalty;
    }


    protected void setLetterToIndex(char[] order) {
        for (int i = 0; i < numEnglishLetters; i++) letterToIndex[i] = -1;
        for (int i = 0; i < numAminoAcids + numWildChars; i++) {
            letterToIndex[order[i] - asciiA] = i;
        }//End of revOrder manufacturing
    }

    public void print(){
        System.out.print("  ");
        for (int orderInd = 0; orderInd < numAminoAcids+numWildChars; orderInd++)
            System.out.print(order[orderInd] + "  ");
        System.out.print("\n");
        for (int row = 0; row < numAminoAcids+numWildChars; row++) {
            for (int col = 0; col < numAminoAcids+numWildChars; col++){
                System.out.printf("%3d", substitutionMatrix[row][col]);
            }
            System.out.print("\n");
        }
    }

    public abstract double minScore();
    public int[][] substitutionMatrix(){
        return substitutionMatrix;
    }

    public int[] letterToIndex(){
        return letterToIndex;
    }

    public char[] order(){
        return order;
    }

    public double gapStartPenalty(){
        return gapStartPenalty;
    }

    public double inGapPenalty(){
        return inGapPenalty;
    }
}