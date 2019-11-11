package meshi.energy.simpleEnergyTerms.repulsion;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParameters;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsParametersList;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.FreeDistance;
import meshi.geometry.FreeDistanceList;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Utils;

import java.util.ArrayList;

/**
 * Flat bottom linear distance constraints function
 */
public class Repulsion extends SimpleEnergyTerm{
    public Repulsion(FreeDistanceList distances,
                                  DistanceConstraintsParametersList minDistance, EnergyInfoElement info){
        super(toArray(distances),minDistance,info);
        if (distances.size() != minDistance.size())
            throw new RuntimeException("This is weird");
        comment = "Repulsion";
        createElementsList(distances,minDistance);
    }
    public EnergyElement createElement(Object baseElement,
                                       Parameters parameters){
        return new RepulsionElement((Distance)baseElement,
                                             ((DistanceConstraintsParameters)parameters).target,weight);
    }

    public boolean evaluatesResidues() {return false;}
    private void  createElementsList(FreeDistanceList distances,DistanceConstraintsParametersList targets){
        elementsList = new ArrayList();
        for (int i = 0; i < distances.size(); i++) {
            Distance distance = distances.get(i);
            DistanceConstraintsParameters target = (DistanceConstraintsParameters) targets.get(i);

            elementsList.add(new RepulsionElement(distance,target.target, weight));
        }
    }
    public void scaleWeight(double factor) {
        super.scaleWeight(factor);
        for (EnergyElement element : elementsList)
            ((RepulsionElement) element).scaleWeight(factor);
    }

}
