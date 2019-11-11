/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

import meshi.sequences.*;

public class SubstitutionScorer implements CellScorer {
	public final static int asciiA=65, numAminoAcids=20, numWildChars=3, numEnglishLetters=26;
	public final static int UP=0, LEFT=1, DIAGONAL=2, MAX=3,numDirections=3, numScores=4;
    private double gapStartPenalty;//added/modified by Tommer, 2.9.14
    private double inGapPenalty;//modified by Tommer, 2.9.14
    private int substitutionMatrix[][]=new int[numAminoAcids+numWildChars][numAminoAcids+numWildChars];//added 1.9.14
    private int letterToIndex[] = new int[numEnglishLetters];//might need deleting 1.9.14

    
    public SubstitutionScorer(SubstitutionMatrix alSchem) {//signature changed Tommer 1.9.14
    	this.gapStartPenalty =alSchem.gapStartPenalty();
        this.inGapPenalty = alSchem.inGapPenalty();
        this.substitutionMatrix=alSchem.substitutionMatrix();//Added
        this.letterToIndex=alSchem.letterToIndex();
    }

    
    public void getScores(Cell cell, Cell[] bestRoutes, double[] scores) {

        DpMatrix matrix = cell.dpMatrix;
        
        if (bestRoutes.length!=numDirections)
        	throw new RuntimeException("cell array needs to be of length 3");
        if (scores.length!=numScores)
        	throw new RuntimeException("scores array needs to be of length 4");
        	
        double internalScore=internalScore(matrix, cell, substitutionMatrix, letterToIndex);     
        
        Cell up = cell.upCell;
        Cell left = cell.leftCell;
        Cell diagonal = cell.diagonalCell;

        double firstScore=Math.max(0, internalScore);
        double[] zeroScores={firstScore,firstScore,firstScore,firstScore};
        if ((up == null) ||(left == null)) {//modified Tommer 9.9.14
            cell.setBack(null, null);
            cell.setBestRoutes(null, null, null);
            for (int iScore = 0; iScore < zeroScores.length; iScore++)
                scores[iScore] = zeroScores[iScore];
                //  scores = zeroScores //modified 23.9.14 Tommer
            return;//Tommer 23.9.14
        }
/*
		if (up == null) {//modified Chen
			cell.setBack(left, left.leftCell);
			cell.setBestRoutes(left.leftCell, left.leftCell, left.leftCell);
			scores[UP] = 0;
			scores[LEFT] = inGapPenalty + left.scoreLeft;
			scores[DIAGONAL] = 0;
			scores[MAX] = scores[LEFT]+internalScore;
			//  scores = zeroScores //modified 23.9.14 Tommer
			return;//chen
		}

		if (left == null) {//modified Chen
			cell.setBack(up, up.upCell);
			cell.setBestRoutes(up.upCell, up.upCell,up.upCell);
			scores[UP] = inGapPenalty + up.scoreUp;
			scores[LEFT] = 0;
			scores[DIAGONAL] = 0;
			scores[MAX] = scores[UP]+internalScore;
			//  scores = zeroScores //modified 23.9.14 Tommer
			return;//chen
		}
	*/
		double leftLeft,leftUp,leftDiagonal;
		if (cell.rowNumber != matrix.sequence1.size()-1)
			leftLeft=inGapPenalty + left.scoreLeft;
		else
			leftLeft=left.scoreLeft;
		leftUp=gapStartPenalty + left.scoreUp;
		leftDiagonal=gapStartPenalty + left.scoreDiagonal;
		double arrayLeft[]={leftLeft, leftUp,leftDiagonal};
		double scoreLeft = maxMultiple(arrayLeft);//maximal score from left

		double upLeft,upUp,upDiagonal;
		if (cell.colNumber != matrix.sequence2.size()-1)
			upUp = inGapPenalty + up.scoreUp;
		else
			upUp = up.scoreUp;
		upLeft=gapStartPenalty + up.scoreLeft;
		upDiagonal=gapStartPenalty + up.scoreDiagonal;
		double arrayUp[]={upLeft, upUp,upDiagonal};
		double scoreUp = maxMultiple(arrayUp);//maximal score from above

		double diagonalLeft,diagonalUp,diagonalDiagonal;
		diagonalUp=diagonal.scoreUp + internalScore;
		diagonalLeft=diagonal.scoreLeft + internalScore;
		diagonalDiagonal=diagonal.scoreDiagonal + internalScore;
		double arrayDiagonal[]={diagonalLeft, diagonalUp,diagonalDiagonal,0,internalScore};
		double scoreDiagonal = maxMultiple(arrayDiagonal);//maximal score from diagonal

		double arrayMax[]={scoreUp, scoreLeft, scoreDiagonal};
        double max = maxMultiple(arrayMax);
        
        
		if (upUp == scoreUp) {
			bestRoutes[UP] = up.upCell;
		} else {
			if (upLeft == scoreUp) {
				bestRoutes[UP] = up.leftCell;
			} else {
				bestRoutes[UP] = up.diagonalCell;
			}
		}
		//System.out.println(bestRoutes[0]);

		if (leftUp == scoreLeft) {
			bestRoutes[LEFT] = left.upCell;
		} else {
			if (leftLeft == scoreLeft) {
				bestRoutes[LEFT] = left.leftCell;
			} else {
				bestRoutes[LEFT] = left.diagonalCell;
			}
		}
		
		if (diagonalUp == scoreDiagonal) {
			bestRoutes[DIAGONAL] = diagonal.upCell;
		} else {
			if (diagonalLeft == scoreDiagonal) {
				bestRoutes[DIAGONAL] = diagonal.leftCell;
			} else {
				if (diagonalDiagonal == scoreDiagonal){
					bestRoutes[DIAGONAL] = diagonal.diagonalCell;
				}

			}
		}
		
		cell.setBestRoutes(bestRoutes[UP], bestRoutes[LEFT], bestRoutes[DIAGONAL]);

        if (scoreUp == max){
        	if (upUp == max)
        		cell.setBack(up,up.upCell);
        	else
        		if (upLeft == max)
        			cell.setBack(up,up.leftCell);
        		else
        			cell.setBack(up,up.diagonalCell);
        }
        else if (scoreLeft == max){
        	if (leftUp == max)
        		cell.setBack(left,left.upCell);
        	else
        		if (leftLeft == max)
        			cell.setBack(left,left.leftCell);
        		else
        			cell.setBack(left,left.diagonalCell);
        }
        else if (scoreDiagonal == max){
        	if (diagonalUp == max)
        		cell.setBack(diagonal,diagonal.upCell);
        	else
        		if (diagonalLeft == max)
        			cell.setBack(diagonal,diagonal.leftCell);
        		else
        			cell.setBack(diagonal,diagonal.diagonalCell);
        }
        
        /*if (cell.rowNumber-cell.colNumber<5&&cell.rowNumber-cell.colNumber>-5&&
        		bestRoutes[0]!=null&&bestRoutes[1]!=null&&bestRoutes[2]!=null){ 
        	
        }*/
        
        scores[UP]=scoreUp;
        scores[LEFT]=scoreLeft;
        scores[DIAGONAL]=scoreDiagonal;
        scores[MAX]=max;//the numbers are odd here...
       // return new double[] {max,scoreUp,scoreLeft,scoreDiagonal}; // new needs to be ELIMINATED!!! Tommer 21.9.14
        
    }
    
    public double maxMultiple(double[] array){
    	try{
    	if (array==null||array.length<1){
    		throw new RuntimeException("null array or empty array"); 
    	}
    	}
    	catch(RuntimeException e){
    		System.out.println("null array or empty array");
    	}
    	
    	double temp=array[0];
    	for (int i=0;i<array.length;i++){
    		temp=Math.max(temp, array[i]);
    	}
    	return temp;
    }
    
    public static double internalScore(DpMatrix matrix, Cell cell, int[][] subMat, int[]revOrder){
    	char rowChar = matrix.rowChar(cell.rowNumber);
        char columnChar = matrix.columnChar(cell.colNumber);
        
        double internalScore;
        
       
        if ((rowChar == SequenceAlignmentCell.GAP_CHAR) |
                (columnChar == SequenceAlignmentCell.GAP_CHAR) |
                (rowChar == SequenceAlignmentCell.WILDCARD_CHAR) |
                (columnChar == SequenceAlignmentCell.WILDCARD_CHAR))
        	internalScore = 0;
        //else internalScore = ((rowChar == columnChar) ? 5 :0);//Edited by Tommer 1.9.14, mismatch penalty was 0
        
        /*the following implementation replaces the above line. it should do the trick
         (of allowing for the use of a 2D matrix score) ASSUMING this is really the only
          relevant place were score is given for matches/mismatches. this is to be checked.  */
        
       else{
        	if(revOrder[(int)rowChar-asciiA]==-1){
        		System.out.println("error: unrecognized residue letter");
        		internalScore =0;
        	}
        	else
        		internalScore =subMat[revOrder[(int)rowChar-65]][revOrder[(int)columnChar-65]];
        		
       } 
        return internalScore;
    }
}

