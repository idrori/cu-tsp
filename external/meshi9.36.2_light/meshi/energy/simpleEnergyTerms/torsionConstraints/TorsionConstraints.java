package meshi.energy.simpleEnergyTerms.torsionConstraints;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.energy.simpleEnergyTerms.distanceMapConstraints.DistanceMapConstraintElement;
import meshi.geometry.*;

import java.util.ArrayList;

public class TorsionConstraints extends SimpleEnergyTerm{
    TorsionList torsions;
    public TorsionConstraints(TorsionList torsions,
                                  TorsionConstraintsParametersList targets, EnergyInfoElement info){
        super(toArray(torsions),targets,info);
        if (torsions.size() != targets.size())
            throw new RuntimeException("This is weird");
        comment = "TorsionConstraints";
        createElementsList(torsions,targets);
        this.torsions = torsions;
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        return null;
    }
    public void scaleWeight(double factor) {
        super.scaleWeight(factor);
        for (EnergyElement element : elementsList) {
            ((TorsionConstraintElement) element).scaleWeight(factor);
        }
    }
    public boolean evaluatesResidues() {return false;}
    public TorsionList torsionList() {return torsions;}

    private void  createElementsList(TorsionList torsions,TorsionConstraintsParametersList targets){
        elementsList = new ArrayList();
        for (int i = 0; i < torsions.size(); i++) {
            Torsion torsion = torsions.get(i);
            TorsionConstraintsParameters target = (TorsionConstraintsParameters) targets.get(i);
            elementsList.add(new TorsionConstraintElement(torsion,target,weight));
        }
    }

    public void setBeta(double beta) {
        for (EnergyElement element : elementsList) {
            ((DistanceMapConstraintElement) element).setBeta(beta);
        }
    }

    public void scaleWeights(double factor) {
        for (EnergyElement element : elementsList) {
            ((TorsionConstraintElement) element).scaleWeight(factor);
        }
    }

    public void report() {
        for (EnergyElement element : elementsList)
            ((TorsionConstraintElement) element).report();
    }
}
