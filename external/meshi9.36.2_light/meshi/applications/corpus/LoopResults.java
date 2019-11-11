/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

public class LoopResults {
    protected double solvEnergy;
    protected double evEnergy;
    protected double propEnergy;
    protected double distEnergy;
    protected double RMS;
    private double[][] bbCoors = null; // for N residues this is a [3N][3] array. Indexing is N,Ca,C,N,Ca,C... and the other index is {x,y,z}
    public int nn, cc, mid;
    public int numSimilar = 0;

    public LoopResults() {
    }

    ;

    public LoopResults(double solvEnergy,
                       double evEnergy,
                       double propEnergy,
                       double distEnergy,
                       double RMS) {
        this.solvEnergy = solvEnergy;
        this.evEnergy = evEnergy;
        this.propEnergy = propEnergy;
        this.distEnergy = distEnergy;
        this.RMS = RMS;
    }

    public void setBBcoors(double[][] coors) {
        bbCoors = coors;
    }

    public double[][] getBBcoors() {
        return bbCoors;
    }

    public double calcRMS(LoopResults other) {
        if (bbCoors[0].length != other.bbCoors[0].length)
            throw new RuntimeException("ERROR: the bbcoors arrays in both instances are not the same length.");
        double totRms = 0.0;
        for (int c = 0; c < bbCoors[0].length; c++)
            totRms +=
                    ((bbCoors[0][c] - other.bbCoors[0][c]) * (bbCoors[0][c] - other.bbCoors[0][c]) +
                            (bbCoors[1][c] - other.bbCoors[1][c]) * (bbCoors[1][c] - other.bbCoors[1][c]) +
                            (bbCoors[2][c] - other.bbCoors[2][c]) * (bbCoors[2][c] - other.bbCoors[2][c]));
        return Math.sqrt(totRms / bbCoors[0].length);
    }
}
