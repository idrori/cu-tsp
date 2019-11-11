/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

public interface CellScorer {
    public void getScores(Cell cell, Cell[] bestRoutes, double[] scores);
}
