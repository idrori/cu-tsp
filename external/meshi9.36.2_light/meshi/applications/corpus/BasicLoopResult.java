/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;


public class BasicLoopResult {
    protected double evEnergy;
    protected double score;
    protected double[][] coors = null; // This pointer refer to the output of the "saveLoopCoordinates()" method
    protected double[][] myPhiPsi = null; // This pointer refer to the output of the "savePhiPsiOfLoop()" method
    protected int[] fragRank = null; // From where in the library did the fragments which makes the solution come?

    public BasicLoopResult() {
    }

    ;

    public BasicLoopResult(double evEnergy, double score, int[] fragRank_p, double[][] coors, double[][] myPhiPsi) {
        this.evEnergy = evEnergy;
        this.score = score;
        this.coors = coors;
        this.myPhiPsi = myPhiPsi;
        fragRank = new int[fragRank_p.length];
        for (int c = 0; c < fragRank_p.length; c++)
            fragRank[c] = fragRank_p[c];
    }

}
