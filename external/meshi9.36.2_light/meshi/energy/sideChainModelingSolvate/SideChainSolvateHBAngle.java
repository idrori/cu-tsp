/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.sideChainModelingSolvate;

import meshi.util.mathTools.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;

/**
 * This class is a truncated form of the corresponding class in the "meshi.energy.solvation" package. It was designed
 * solely for accelerating the application SCMOD (concurrent sidechain modeling), and should not be used for other
 * purposes.
 * <p/>
 * Please do not use it!
 */


public class SideChainSolvateHBAngle {

    private DistanceMatrix dm;
    private final double sigmoidBeginsWithH;
    private final double sigmoidEndsWithH;
    private final double sigmoidBeginsNoH;
    private final double sigmoidEndsNoH;
    private double sigmoidBegins;
    private double sigmoidEnds;
    private double hbAngScore;    // This is a variable solvation uses, and has a getter method
    private double sigmCosAng1;
    private double sigmCosAng2;
    private final Sigma sigma = new Sigma();
    /**
     * The constructor parameters:
     * dm - The distance matrix used in this object. All the getAtoms given as parameters to the 'updateAndEvaluateAtoms1234'
     * method must be represented in this object.
     * sigmoidBeginsWithH,sigmoidEndsWithH - The transition angles (in degrees) for the HB sigmoid when the hydrogen
     * atom is known explicitly. Above sigmoidEndsWithH the sigmoid is given a value of 1.0 . Bellow sigmoidBeginsWithH
     * the sigmoid is given a value of 0.0 . In between it raises smoothly by cubic spline.
     * sigmoidBeginsNoH,sigmoidEndsNoH - The same as above, only for HB angle sigmoids when the hydrogen atom is not
     * given explicitly.
     */
    public SideChainSolvateHBAngle(DistanceMatrix dm,
                                   double sigmoidBeginsWithH,
                                   double sigmoidEndsWithH,
                                   double sigmoidBeginsNoH,
                                   double sigmoidEndsNoH) {
        this.dm = dm;
        this.sigmoidBeginsNoH = Math.cos(sigmoidEndsNoH * Math.PI / 180.0); // note the switch end-begin. cos is a monotonicly DECREASING function
        this.sigmoidEndsNoH = Math.cos(sigmoidBeginsNoH * Math.PI / 180.0);
        this.sigmoidBeginsWithH = Math.cos(sigmoidEndsWithH * Math.PI / 180.0);
        this.sigmoidEndsWithH = Math.cos(sigmoidBeginsWithH * Math.PI / 180.0);
    }

    public final double hbAngScore() {
        return hbAngScore;
    }


    public final void updateAndEvaluateAtoms1234(Atom a1, Atom a2, Atom a3, Atom a4) {
        if (a1.nowhere() | a2.nowhere() | a3.nowhere() | a4.nowhere()) return;
        if (a2.type().isHydrogen() || a3.type().isHydrogen()) {
            sigmoidBegins = sigmoidBeginsWithH;
            sigmoidEnds = sigmoidEndsWithH;
        } else {
            sigmoidBegins = sigmoidBeginsNoH;
            sigmoidEnds = sigmoidEndsNoH;
        }

        updateAndEvaluateAtoms123(a1, a2, a3);
        updateAndEvaluateAtoms234(a2, a3, a4);

        hbAngScore = sigmCosAng1 * sigmCosAng2;
    }

    private final void updateAndEvaluateAtoms123(Atom a1, Atom a2, Atom a3) {
        Distance dis1 = dm.distance(a1, a2); // distance isOn pair (a1-a2)
        if (dis1 == null) dis1 = new DistanceMirror(dm.distance(a2, a1));
        Distance dis2 = dm.distance(a3, a2); // distance isOn pair (a3-a2)
        if (dis2 == null){
            if (dm.distance(a2, a3)== null) {
                throw new RuntimeException("This is weird "+a2.distanceFrom(a3)+"\n"+a2+"\n"+a3);
            }
            dis2 = new DistanceMirror(dm.distance(a2, a3));
        }
        double cosAng = dis1.dDistanceDx() * dis2.dDistanceDx() + dis1.dDistanceDy() * dis2.dDistanceDy() + dis1.dDistanceDz() * dis2.dDistanceDz();

        sigma.sigma(1.0 + cosAng, 2, 1.0 + sigmoidBegins, 1.0 + sigmoidEnds, 1.0, 0.0);
        sigmCosAng1 = sigma.s();
    }

    private final void updateAndEvaluateAtoms234(Atom a2, Atom a3, Atom a4) {
        Distance dis1 = dm.distance(a2, a3); // distance isOn pair (a2-a3)
        if (dis1 == null) dis1 = new DistanceMirror(dm.distance(a3, a2));
        Distance dis2 = dm.distance(a4, a3); // distance isOn pair (a4-a3)
        if (dis2 == null) dis2 = new DistanceMirror(dm.distance(a3, a4));

        double cosAng = dis1.dDistanceDx() * dis2.dDistanceDx() + dis1.dDistanceDy() * dis2.dDistanceDy() + dis1.dDistanceDz() * dis2.dDistanceDz();

        sigma.sigma(1.0 + cosAng, 2, 1.0 + sigmoidBegins, 1.0 + sigmoidEnds, 1.0, 0.0);
        sigmCosAng2 = sigma.s();
    }
}	
