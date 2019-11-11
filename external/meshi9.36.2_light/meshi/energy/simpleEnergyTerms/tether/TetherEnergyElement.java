/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.tether;

import meshi.energy.EnergyElement;
import meshi.geometry.Coordinates;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;


public class TetherEnergyElement extends EnergyElement {
    protected Atom atom, evaluatedAtom;
    protected double force, forceSave;
    protected double evaluatedX;
    protected double evaluatedY;
    protected double evaluatedZ;
    //protected boolean frozen;
    protected FreeDistance distance;
    private double reliability;

    public double reliability() {
        return reliability;
    }

    public void setReliability(double reliability) {
//        if (reliability != 0) throw new RuntimeException("This is weird");
        this.reliability = reliability;
    }

    public TetherEnergyElement() {
    }

    public TetherEnergyElement(Atom inputAtom, double x, double y, double z, double weight) {
        evaluatedX = x;
        evaluatedY = y;
        evaluatedZ = z;

        atom = inputAtom;
        atoms = new AtomList(atom.molecularSystem);
        atoms.add(atom);
        reset(evaluatedX, evaluatedY, evaluatedZ);
        force = weight;
        forceSave = weight;
        frozen = atom.frozen();
        reliability = 1;
    }


    protected void reset() {
        reset(atom.x(), atom.y(), atom.z());
    }

    protected void reset(double x, double y, double z) {
        MolecularSystem copyMS = new MolecularSystem();
        evaluatedAtom = new Atom("eval", null, AtomType.XXX, new Coordinates(x, y, z), 1, 0.0,copyMS);
        distance = new FreeDistance(atom, evaluatedAtom);
    }

    protected void setAtoms() {
    }


    public void remove() {
        force = 0;
    }

    public void restoreWeight() {
        force = forceSave;
    }

    public double evaluate() {
        double d;
        double deDd;
        double deDx;
        double deDy;
        double deDz;
        double energy;
        double dMinusEpsilon;

        if (frozen) return 0;
        distance.update();
        d = distance.distance();
        if (d == 0) return 0.0;

        energy = d * d * force * reliability;
        deDd = d * force * 2 * reliability;


        deDx = deDd * distance.dDistanceDx();
        deDy = deDd * distance.dDistanceDy();
        deDz = deDd * distance.dDistanceDz();
        if (!atom.frozen()) {
            atom.addToFx(-1 * deDx); // force = -derivative
            atom.addToFy(-1 * deDy);
            atom.addToFz(-1 * deDz);
        }

        return energy;
    }


    /*
private double distance(Atom atom, double X,double Y, double Z){

return Math.sqrt((atom.x()-X)*(atom.x()-X)+(atom.y()-Y)*(atom.y()-Y)+(atom.z()-Z)*(atom.z()-Z));

}                                             */

    public void update() {
    }

    public String toString() {
        return atom.toString();
    }


}
