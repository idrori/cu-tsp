package meshi.energy.cooperativeDistanceMapConstraints;


import meshi.energy.EnergyElement;
import meshi.geometry.FreeDistance;
import meshi.geometry.FreeDistanceList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;

import java.util.ArrayList;

/**
 * A flat bottom linear distance constraint element
 */
public class CooperativeDistanceMapConstraintElement extends EnergyElement{
    public final double EPSILON = 0.1;
    private ArrayList<double[]> targetsAndRanges;
    private FreeDistanceList distances;

    private double weight;
    private double beta;
    double[] betas;
    double[] diffs;
    public CooperativeDistanceMapConstraintElement(FreeDistanceList distances, ArrayList<double[]> targetsAndRanges, double beta, double weight){
        this.distances = distances;
        this.targetsAndRanges = targetsAndRanges;
        this.weight = weight;
        this.beta = beta;
        betas = new double[distances.size()];
        diffs = new double[distances.size()];
        setAtoms();
    }
    protected void setAtoms() {
        atoms = new AtomList(distances.get(0).atom1());
        for (FreeDistance distance : distances) {
            Atom atom1 = distance.atom1();
            Atom atom2 = distance.atom2();
            if (! atoms().contains(atom1)) atoms.add(atom1);
            if (! atoms().contains(atom2)) atoms.add(atom2);
        }

    }

    public double evaluate() {
        double energy = -weight;
        double factor;
        for (int i = 0; i < distances.size(); i++) {
            FreeDistance distance = distances.get(i);
            double       dis      = distance.distance();
            double       target   = targetsAndRanges.get(i)[0];
            double       range    = targetsAndRanges.get(i)[1];
            double diff;
            if (dis >= target + range/2)
                diff = dis - (target + range/2);
            else {
                if (dis <= target - range/2)
                    diff = dis - (target - range/2);
                else {
                    diff = 0;
                }
            }
            if (diff > 0)
                factor = 1;
            else factor = 4/beta;
            double       exp      = Math.exp(-beta*factor*diff*diff);
            energy = energy * exp;
            diffs[i] = diff;
            betas[i]  = beta*factor;
        }
        for (int i = 0; i < distances.size(); i++) {
            FreeDistance distance = distances.get(i);
            double diff = diffs[i];
            double beta  = betas[i];
            double dEdDistance = -energy * 2 * beta * diff;
            Atom atom1 = distance.atom1();
            Atom atom2 = distance.atom2();
            atom1.addToFx(-dEdDistance * distance.dDistanceDx());
            atom1.addToFy(-dEdDistance * distance.dDistanceDy());
            atom1.addToFz(-dEdDistance * distance.dDistanceDz());
            atom2.addToFx(dEdDistance * distance.dDistanceDx());
            atom2.addToFy(dEdDistance * distance.dDistanceDy());
            atom2.addToFz(dEdDistance * distance.dDistanceDz());
        }
        return energy;
    }
}
