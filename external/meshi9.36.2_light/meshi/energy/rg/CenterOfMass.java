/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.Atom;
import meshi.util.*;

/**
 * Center of mass.
 */
public class CenterOfMass implements Updateable {
    private double x;
    private double y;
    private double z;
    /*
     * Used to avoid reupdate in the same minimization step
     */
    private int numberOfUpdates = 0;


    public final double getX() {
        return x;
    }

    public final double getY() {
        return y;
    }

    public final double getZ() {
        return z;
    }

    private AtomList atomList;
    private int size;
    private double invSize;

    public final double invSize() {
        return invSize;
    }

    public final int size() {
        return size;
    }

    public CenterOfMass(AtomList atomList) {
        this.atomList = atomList;
        size = atomList.size();
        invSize = 1.00 / size;
        update();
    }

    public void update(int numberOfUpdates) {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            this.numberOfUpdates++;
            update();
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with HbondList.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }

    private void update() {
        size = atomList.size();
        x = 0;
        y = 0;
        z = 0;
        if (size == 0) return;
        for (Atom atom : atomList) {
            x += atom.x();
            y += atom.y();
            z += atom.z();
        }
        x /= size;
        y /= size;
        z /= size;
    }

    public Atom atomAt(int i) {
        return atomList.get(i);
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }
}

