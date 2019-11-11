package meshi.sequences.aligner;

import meshi.util.file.MeshiLineReader;

import java.io.IOException;

public class IdentityMatrix extends SubstitutionMatrix{
    public final static double MIN_SCORE = -10.5;
	public final static int asciiA=65;
    public static final double GAP_START_PENALTY = -5.0;
    public static final double IN_START_PENALTY = -0.5;
//    public static final double GAP_START_PENALTY = 0.0;
//    public static final double IN_START_PENALTY = 0.0;


    public double minScore() {
        return MIN_SCORE;
    }
    public IdentityMatrix()  {
        super(GAP_START_PENALTY,IN_START_PENALTY);
        for (int iRow = 0; iRow < substitutionMatrix.length;iRow++)
            for (int jColumn = 0; jColumn < substitutionMatrix.length;jColumn++)
                if (iRow == jColumn) substitutionMatrix[iRow][jColumn]=1;
                else substitutionMatrix[iRow][jColumn]=-12;
        String temp = "ARNDCQEGHILKMFPSTWYVBZX";
        order = temp.toCharArray();
        setLetterToIndex(order);
    }
}
