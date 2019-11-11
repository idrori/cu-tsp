package meshi.energy.simpleEnergyTerms.repulsion;


import meshi.energy.EnergyElement;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParameters;
import meshi.geometry.Distance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.Utils;

/**
 * A flat bottom linear distance constraint element
 */
public class RepulsionElement extends EnergyElement{
    public final double EPSILON = 0.1;
    private double minD;
    private Distance distance;
    private Atom atom1,atom2;
    private double weight;
    private boolean debug;

    public RepulsionElement(Distance distance, double minD, double weight){
        this.minD = minD;
        this.distance = distance;
        this.weight = weight;
        atom1 = distance.atom1();
        atom2 = distance.atom2();
        setAtoms();
        if (((atom1.residueNumber() == 61) & atom1.name().equals("O")) & ((atom2.residueNumber() == 62) & atom2.name().equals("CB")))
            debug = true;
        else debug = false;
    }
    protected void setAtoms() {
        atoms = new AtomList(2,distance.atom1().molecularSystem);
        atoms.add(atom1);
        atoms.add(atom2);
        atom1 = distance.atom1();
        atom2 = distance.atom2();
    }

    public double evaluate() {
        double dis = distance.distance();
//        if (debug)
//            Utils.printDebug(this, "xxxxxxxxxxxxxxxxxxxxxxxxxxx "+dis+" "+minD);
        double dx  = distance.dDistanceDx();
        double dy  = distance.dDistanceDy();
        double dz  = distance.dDistanceDz();
        double disMinusMinD;
        if (dis > minD)
            return 0;
        disMinusMinD = dis - minD;
//        if (debug)
//            Utils.printDebug(this, "yyyyyyyyyy "+dis+" "+minD);
        double disMinusMinD2 = disMinusMinD * disMinusMinD;
        double energy = weight * disMinusMinD2;
        double force = -2 * weight * disMinusMinD;
            double dForceDx = force * dx;
            double dForceDy = force * dy;
            double dForceDz = force * dz;
            atom1.addToFx(dForceDx);
            atom1.addToFy(dForceDy);
            atom1.addToFz(dForceDz);
            atom2.addToFx(-dForceDx);
            atom2.addToFy(-dForceDy);
            atom2.addToFz(-dForceDz);
        return energy;
    }
    public void scaleWeight(double factor) {
        weight *= factor;
    }


}
