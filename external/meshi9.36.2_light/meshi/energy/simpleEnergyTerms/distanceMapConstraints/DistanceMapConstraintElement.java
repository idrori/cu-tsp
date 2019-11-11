package meshi.energy.simpleEnergyTerms.distanceMapConstraints;


import meshi.energy.EnergyElement;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParameters;
import meshi.geometry.Distance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.Utils;

/**
 * A flat bottom linear distance constraint element
 */
public class DistanceMapConstraintElement extends EnergyElement{
    public final double EPSILON = 0.1;
    private double targetBottom, targetTop;
    private Distance distance;
    private Atom atom1,atom2;
    private double weight;
    private double beta;
    double buffer = 0;
    int separation;
    double scale;
    boolean on;

    public DistanceMapConstraintElement(Distance distance, DistanceConstraintsParameters targetAndRange, double weight){
        this.distance = distance;
        double target = targetAndRange.target;
        double range = targetAndRange.range;
        this.targetBottom = target - range/2;
        this.targetTop = target + range/2;
        this.weight = weight;
        atom1 = distance.atom1();
        atom2 = distance.atom2();
        separation = atom2.residueNumber() - atom1.residueNumber();
        beta = 1;
        scale = calcScale(beta);
//        if (scale < 0.001)
//            on = true;
//        else
//            on = false;
//        if (target < 0) {
//            on = false;
//            targetBottom = -target - range/2;
//            targetTop    = -target + range/2;
//        }
        on = true;
        setAtoms();
    }
    private double calcScale(double beta) {
        return Math.exp(-beta*(separation - 1));
    }
    public void setBeta(double beta) {
        this.beta = beta;
        scale = calcScale(beta);
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
        double disMinusTarget;
        if (dis >= targetTop)
            disMinusTarget = dis - targetTop;
        else {
            if (dis <= targetBottom)
                disMinusTarget = dis - targetBottom;
            else
                return 0;
        }
        if (Math.abs(disMinusTarget) < 0.4)
            DistanceMapConstraints.distanceCount++;

        double disMinusTarget2;
        double disMinusTarget3;
        double disMinusTarget4;
        double disMinusBottom3PlusEpsilon;
        double energy;
        double force;

        if (disMinusTarget <= 0) {
            disMinusTarget2 = disMinusTarget*disMinusTarget;
            energy = 1000*weight * disMinusTarget2 ;
            force = -1000*weight * 2 * disMinusTarget ;
        }
        else {
            disMinusTarget = disMinusTarget - buffer;
            if (on & (disMinusTarget > 0)) {
//            if (disMinusTarget > 0) {
                disMinusTarget2 = disMinusTarget*disMinusTarget;
                disMinusTarget3 = disMinusTarget2*disMinusTarget;
                disMinusTarget4 = disMinusTarget2*disMinusTarget2;
                disMinusBottom3PlusEpsilon = EPSILON +disMinusTarget3;
                energy = scale * weight * disMinusTarget4 / disMinusBottom3PlusEpsilon;
                force = scale * (-weight * disMinusTarget3 / disMinusBottom3PlusEpsilon) *
                        (4 - 3 * disMinusTarget3 / disMinusBottom3PlusEpsilon);
            } else {
//                if (disMinusTarget <= 0)
//                    on = false;
                energy = 0;
                force = 0;
            }
        }
//            if (target <= 8) {
//                energy = energy - 0.1 * scale * weight * Math.exp(-0.1*disMinusTarget2) + 0.1 * scale * weight;
//                force  = force  + 0.1 * scale * weight * Math.exp(-0.1*disMinusTarget2) * -0.2 * disMinusTarget ;
//            }
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


}
