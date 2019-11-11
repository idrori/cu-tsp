package meshi.energy.simpleEnergyTerms.distanceMapConstraints;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintElement;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParameters;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParametersList;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.energy.simpleEnergyTerms.torsionConstraints.TorsionConstraintElement;
import meshi.geometry.Distance;
import meshi.geometry.FreeDistance;
import meshi.geometry.FreeDistanceList;
import meshi.util.Utils;

import java.util.ArrayList;

/**
 * Flat bottom linear distance constraints function
 */
public class DistanceMapConstraints extends SimpleEnergyTerm{
    protected static int distanceCount;
    private int numberOfDistance;
    private int numberOfUpdates;
    public DistanceMapConstraints(FreeDistanceList distances,
                                  DistanceConstraintsParametersList targets, EnergyInfoElement info){
        super(toArray(distances),targets,info);
        if (distances.size() != targets.size())
            throw new RuntimeException("This is weird");
        comment = "DistanceMapConstraints";
        createElementsList(distances,targets);
        numberOfDistance = 0;
        numberOfUpdates = 0;
    }
    public EnergyElement createElement(Object baseElement,
                                       Parameters parameters){
        double weight = 1;
        return new DistanceMapConstraintElement((Distance)baseElement,
                                             ((DistanceConstraintsParameters)parameters),weight);
    }

    public boolean evaluatesResidues() {return false;}
    @Override
    public EnergyInfoElement evaluate() throws EvaluationException {
            distanceCount = 0;
            EnergyInfoElement energyInfoElement = super.evaluate();
            if (distanceCount != numberOfDistance) {
                numberOfDistance = distanceCount;
                numberOfUpdates++;
                if (numberOfUpdates % 100 == 0)
                    Utils.println("# of distances: "+ numberOfDistance);
            }
            return energyInfoElement;
    }
    private void  createElementsList(FreeDistanceList distances,DistanceConstraintsParametersList targets){
        elementsList = new ArrayList();
        for (int i = 0; i < distances.size(); i++) {
            Distance distance = distances.get(i);
            DistanceConstraintsParameters target = (DistanceConstraintsParameters) targets.get(i);
            elementsList.add(new DistanceMapConstraintElement(distance,target, weight));
        }
    }



    public void setBeta(double beta) {
        for (EnergyElement element : elementsList) {
            ((DistanceMapConstraintElement) element).setBeta(beta);
        }
    }
}
