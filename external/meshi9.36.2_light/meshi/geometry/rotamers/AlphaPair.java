package meshi.geometry.rotamers;

import meshi.geometry.AlphaCoordinates;

public class AlphaPair {
    private AlphaCoordinates alpha0, alpha1;
    public AlphaPair(AlphaCoordinates alpha0, AlphaCoordinates alpha1) {
        if (alpha0.isValid())
            throw new RuntimeException("Alpha 0 is invalid.");
        if (alpha1.isValid())
            throw new RuntimeException("Alpha 1 is invalid.");
        this.alpha0 = alpha0;
        this.alpha1 = alpha1;
    }

    public double getAngle0() {return alpha0.getAngle0();}
    public double getAngle1() {return alpha0.getAngle1();}
    public double getAngle2() {return alpha1.getAngle1();}
    public double getAlpha0() {return alpha0.getAlpha();}
    public double getAlpha1() {return alpha1.getAlpha();}
}
