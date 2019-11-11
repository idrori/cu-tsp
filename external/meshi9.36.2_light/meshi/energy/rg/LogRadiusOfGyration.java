/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.Updateable;

/**
 *
 */
public class LogRadiusOfGyration implements Updateable {
    private double[] dLogRGdX, dLogRGdY, dLogRGdZ;

    public final double dLogRGdCMx() {
        return rg.dRGdCMX() / rg.radiusOfGyration();
    }

    public final double dLogRGdCMy() {
        return rg.dRGdCMy() / rg.radiusOfGyration();
    }

    public final double dLogRGdCMz() {
        return rg.dRGdCMz() / rg.radiusOfGyration();
    }

    private double logRG;

    public double logRG() {
        return logRG;
    }

    private int size;

    public final double dLogRGdX(int index) {
        return dLogRGdX[index];
    }

    public final double dLogRGdY(int index) {
        return dLogRGdY[index];
    }

    public final double dLogRGdZ(int index) {
        return dLogRGdZ[index];
    }

    public final double dLogRGdX(Atom atom) {
        for (int i = 0; i < size; i++)
            if (atomAt(i) == atom) return dLogRGdX(i);
        return -99999;
    }
    //
    private double radiusOfGyration, dRadiusOfGyration;
    //
    /*
     * Used to avoid reupdate in the same minimization step
     */
    private int numberOfUpdates = 0;
    //
    public final RadiusOfGyration rg;


    public LogRadiusOfGyration(RadiusOfGyration rg) {
        this.rg = rg;
        this.size = rg.size();
        //
        dLogRGdX = new double[size];
        dLogRGdY = new double[size];
        dLogRGdZ = new double[size];
    }

    public void update(int numberOfUpdates) {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            this.numberOfUpdates++;
            rg.update(numberOfUpdates);
            update();
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with HbondList.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }

    private void update() {
        radiusOfGyration = rg.radiusOfGyration();
        logRG = Math.log(radiusOfGyration);
        double invRadiusOfGyration = 1 / radiusOfGyration;
        for (int i = 0; i < size; i++) {
            dLogRGdX[i] = invRadiusOfGyration * rg.dRGdX(i);
            dLogRGdY[i] = invRadiusOfGyration * rg.dRGdY(i);
            dLogRGdZ[i] = invRadiusOfGyration * rg.dRGdZ(i);

        }
    }

    public final AtomList atoms() {
        return rg.atoms();
    }

    public final int size() {
        return rg.size();
    }

    public final Atom atomAt(int index) {
        return rg.atomAt(index);
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
        rg.setNumberOfUpdates(numberOfUpdates);
    }
}
