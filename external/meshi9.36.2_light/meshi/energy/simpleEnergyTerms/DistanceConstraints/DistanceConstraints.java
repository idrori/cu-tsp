package meshi.energy.simpleEnergyTerms.DistanceConstraints;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.FreeDistance;
import meshi.geometry.FreeDistanceList;

import java.util.ArrayList;

/**
 * Flat bottom linear distance constraints function
 */
public class DistanceConstraints extends SimpleEnergyTerm{
    private FreeDistance[] distances;
    public DistanceConstraints(FreeDistanceList distances,
                               DistanceConstraintsParametersList targets,EnergyInfoElement info){
        super(toArray(distances),targets,info);
        if (distances.size() != targets.size())
            throw new RuntimeException("This is weird");
        comment = "DistanceConstraints";
        createElementsList(distances,targets);
    }
    public EnergyElement createElement(Object baseElement,
                                       Parameters parameters){
        return new DistanceConstraintElement((Distance)baseElement,
                                             ((DistanceConstraintsParameters)parameters).target,weight);
    }

    public boolean evaluatesResidues() {return false;}
    private void  createElementsList(FreeDistanceList distances,DistanceConstraintsParametersList targets){
        elementsList = new ArrayList();
        for (int i = 0; i < distances.size(); i++) {
            Distance distance = distances.get(i);
            DistanceConstraintsParameters target = (DistanceConstraintsParameters) targets.get(i);
            elementsList.add(new DistanceConstraintElement(distance,target.target,weight));
        }
    }
}
