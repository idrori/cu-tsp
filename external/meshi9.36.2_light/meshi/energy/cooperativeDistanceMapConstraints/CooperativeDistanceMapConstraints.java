package meshi.energy.cooperativeDistanceMapConstraints;

import meshi.energy.*;
import meshi.geometry.FreeDistance;
import meshi.geometry.FreeDistanceList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Utils;

import java.util.ArrayList;

/**
 * Flat bottom linear distance constraints function
 */
public class CooperativeDistanceMapConstraints extends AbstractEnergy{

    protected double beta = 0.001  ;
    FreeDistanceList freeDistances;
    String mode = null;
    ArrayList<CooperativeDistanceMapConstraintElement> elements;
    public CooperativeDistanceMapConstraints(Protein protein, FreeDistanceList allDistances,
                                             ArrayList<double[]> targets, EnergyInfoElement info, String mode){
        super(toArray(allDistances),info, EnergyType.DIFFERENTIAL);
        this.mode = mode;
        if (allDistances.size() != targets.size())
            throw new RuntimeException("This is weird");
        comment = "CoopDistanceConstraints";
        createElementsList(protein,allDistances, targets);
        this.freeDistances = allDistances;
    }

    private void  createElementsList(Protein protein, FreeDistanceList distances,ArrayList<double[]> targets){
        elements = new ArrayList<>();
        for (Residue residue : protein.chain()) {
            FreeDistanceList newDistanceList = new FreeDistanceList();
            ArrayList<double[]> newTargets = new ArrayList<>();
            if (!residue.dummy() && (!residue.ca().nowhere())) {
                Atom atom;
                if (mode.equals("ca") | (residue.cb() == null))
                    atom = residue.ca();
                else if (mode.equals("cb"))
                        atom = residue.cb();
                else throw new RuntimeException("This is weird");

                for (int i = 0; i < distances.size(); i++) {
                    FreeDistance distance = distances.get(i);
                    if ((distance.atom1() == atom) | (distance.atom2() == atom)) {
                        double dis = distance.distance();
                        double[] target = targets.get(i);
                        if ((dis < target[0]) | (calc(dis, target[0], beta) > calc(0,15, beta))){
                            newDistanceList.add(distance);
                            newTargets.add(target);
                        }
                    }
                }
                if (newDistanceList.size() > 0)
                    elements.add(new CooperativeDistanceMapConstraintElement(newDistanceList, newTargets, beta, weight));
            }
        }
    }

    public EnergyInfoElement evaluate(){
        double energy = 0;
//        for (FreeDistance distance : freeDistances)
//            distance.update();
        for (CooperativeDistanceMapConstraintElement element : elements) {
            energy = energy + element.evaluate();
        }
        info.setValue(energy);
        return info;
    }

    public static double calc(double x, double target, double beta) {
        double d = x - target;
        return Math.exp(-beta*d*d);
    }
    public void test(TotalEnergy energ, Atom atom) {}
    public void evaluateAtoms() {}
    public boolean evaluatesResidues() {return false;}



}
