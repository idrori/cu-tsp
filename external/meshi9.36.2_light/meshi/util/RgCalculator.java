/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.molecularElements.atoms.*;

import java.util.*;

public class RgCalculator implements Updateable {
    private int numberOfUpdates = 0;
    private AtomList atomList;
    private double cmx, cmy, cmz, rg;

    public RgCalculator(AtomList atomList) {
        this.atomList = atomList;
        update();
    }

    public void update(int numberOfUpdates) throws UpdateableException {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            this.numberOfUpdates++;
            update();
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with RgCalculator.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }

    public void update() {
        cmx = cmy = cmz = rg = 0;
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            cmx += atom.x();
            cmy += atom.y();
            cmz += atom.z();
        }
        double size = atomList.size();
        cmx /= size;
        cmy /= size;
        cmz /= size;
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            double dx = cmx - atom.x();
            double dy = cmy - atom.y();
            double dz = cmz - atom.z();
            rg += dx * dx + dy * dy + dz * dz;
        }
        rg /= size;
        rg = Math.sqrt(rg);
    }

    public double rg() {
        return rg;
    }

    public double cmx() {
        return cmx;
    }

    public double cmy() {
        return cmy;
    }

    public double cmz() {
        return cmz;
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }
}
 
