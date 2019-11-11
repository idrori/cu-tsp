/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

import meshi.sequences.*;

public class StructureBasedMatrix implements CellScorer {
    public final double gapPenalty;
    public final double[][] matrixScores;

    public double getScoreUp(Cell cell){//Added by Tommer 2.9.14 to fit with new CellScorer
    	return 0.;
    }
    public double getScoreLeft(Cell cell){
    	return 0.;
    }//end of addition.
    
    public StructureBasedMatrix(double gapPenalty, double[][] scores) {
        this.gapPenalty = gapPenalty;
        this.matrixScores = scores;
    }

    public void getScores(Cell cell, Cell[] bestRoutes, double[] scores) {//modified by Tommer 28.9.14, needs MORE consideration
        DpMatrix matrix = cell.dpMatrix;

        int row = cell.rowNumber;
        int column = cell.colNumber;
        double score = matrixScores[row][column];//modified Tommer 28.9.14

        Cell up = cell.upCell;
        Cell left = cell.leftCell;
        Cell upLeft = cell.diagonalCell;

        if ((up == null) & (left == null)) {
            cell.setBack(null, null);
            for (int i=0;i<4;i++)scores[i]=0;//modified 28.9.14 Tommer
            return;
        }
        if (up == null) {
        	if(left.leftCell == null){
        		cell.setBack(left, null);
        	}
        	else{
        		cell.setBack(left, left);
        	}            
        	for (int i=0;i<4;i++)scores[i]=0;//modified 28.9.14 Tommer
            return;
        }
		if (left == null) {
			if (up.upCell == null) {
				cell.setBack(up, null);
			} else {
				cell.setBack(up, up);
			}
			for (int i=0;i<4;i++)scores[i]=0;//modified 28.9.14 Tommer
            return;
		}

        double scoreUp = (matrix.sequence2.size() == column) ? up.scoreDiagonal : gapPenalty + up.scoreDiagonal;
        double scoreLeft = (matrix.sequence1.size() == row) ? left.scoreDiagonal : gapPenalty + left.scoreDiagonal;
        double scoreUpLeft = upLeft.scoreDiagonal + score;

        double max = Math.max(scoreUp, scoreLeft);
        max = Math.max(max, scoreUpLeft);

        if (scoreUp == max) cell.setBack(up,null);//this is not right! TEMPORARY SOLUTION! Tommer 4.9.14
        else if (scoreLeft == max) cell.setBack(left,null);
        else if (scoreUpLeft == max) cell.setBack(upLeft,null);
        
        scores[0]=max;
        for (int i=1;i<4;i++)scores[i]=0;//modified 28.9.14 Tommer
        return;
    }
}
