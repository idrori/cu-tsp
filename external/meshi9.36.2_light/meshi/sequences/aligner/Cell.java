/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 29/07/2005
 * Time: 15:39:55
 * To change this template use File | Settings | File Templates.
 */
public class Cell {
    public final int rowNumber;
    public final int colNumber;
    public final static int UP=0, LEFT=1, DIAGONAL=2, MAX=3;
    public final Cell upCell, diagonalCell, leftCell;//the next cell will be either the cell in the right of this cell or the cell in bottom of this one
    protected Cell back, nextBack;//nextBack added by Tommer 4.9.14
    public final double maxScore, scoreUp, scoreLeft, scoreDiagonal;
    public final DpMatrix dpMatrix;
    protected Cell[] bestRoutes;
    private double[] scores;
    Cell(int i, int j, DpMatrix mat) {
        rowNumber = i;
        colNumber = j;
        dpMatrix = mat;
        bestRoutes=mat.bestRoutesMatrix()[i][j];
        scores=mat.scoresMatrix()[i][j];
        //the initation of the neighbour cells
        if (rowNumber == 0) {
            upCell = null;
        } else {
            upCell = dpMatrix.getCell(rowNumber - 1, colNumber);
        }

        if (rowNumber == 0 || colNumber == 0) {
            diagonalCell = null;
        } else {
            diagonalCell = dpMatrix.getCell(rowNumber - 1, colNumber - 1);
        }

        if (colNumber == 0) {
            leftCell = null;
        } else {
            leftCell = dpMatrix.getCell(rowNumber, colNumber - 1);
        }
        //the score for the cell
        dpMatrix.cellScorer.getScores(this, bestRoutes, scores);
        scoreUp =scores[UP];
        scoreLeft = scores[LEFT];
        scoreDiagonal = scores[DIAGONAL];//modified Tommer 4.9.14
        maxScore=scores[MAX];
    }

    public void setBack(Cell back, Cell nextBack) {//modified Tommer 4.9.14
        this.back = back;
        this.nextBack=nextBack;
    }
    
    public void setBestRoutes(Cell upNextBack, Cell leftNextBack, Cell diagonalNextBack){//created Tommer 4.9.14
    	bestRoutes[UP]=upNextBack;
    	bestRoutes[LEFT]=leftNextBack;
    	bestRoutes[DIAGONAL]=diagonalNextBack;
    }

    public Cell getBack() {
        return back;
    }

    public String toString() {
        return "cell " + rowNumber + " " + colNumber + " " + maxScore;
    }

}
