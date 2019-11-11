package meshi.energy.simpleEnergyTerms.DistanceConstraints;


import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;

/**
 * A flat bottom linear distance sonstrint element
 */
public class DistanceConstraintElement extends EnergyElement{
    public final double EPSILON = 0.1;
    private double bottom;
    private Distance distance;
    private Atom atom1,atom2;
    private double weight;
    public DistanceConstraintElement(Distance distance, double bottom, double weight){
        this.distance = distance;
        this.bottom = bottom;
        this.weight = weight;
        atom1 = distance.atom1();
        atom2 = distance.atom2();
        setAtoms();
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
        double dx  = distance.dDistanceDx();
        double dy  = distance.dDistanceDy();
        double dz  = distance.dDistanceDz();
        if (dis <= bottom) return 0;

        double disMinusBottom = dis-bottom;
        double disMinusBottom2 = disMinusBottom*disMinusBottom;
        double disMinusBottom3 = disMinusBottom2*disMinusBottom;
        double disMinusBottom4 = disMinusBottom2*disMinusBottom2;
        double disMinusBottom3PlusEpsilon = EPSILON+disMinusBottom3;
        double disMinusBottom3PlusEpsilon2 = disMinusBottom3PlusEpsilon*
                                             disMinusBottom3PlusEpsilon;
        double energy = weight*disMinusBottom4/disMinusBottom3PlusEpsilon;
        double force  = (-weight*disMinusBottom3/disMinusBottom3PlusEpsilon)*
                        (4 - 3*disMinusBottom3/disMinusBottom3PlusEpsilon);
        double dForceDx = force*dx;
        double dForceDy = force*dy;
        double dForceDz = force*dz;
        atom1.addToFx(dForceDx);
        atom1.addToFy(dForceDy);
        atom1.addToFz(dForceDz);
        atom2.addToFx(-dForceDx);
        atom2.addToFy(-dForceDy);
        atom2.addToFz(-dForceDz);
        return energy;
    }
}
