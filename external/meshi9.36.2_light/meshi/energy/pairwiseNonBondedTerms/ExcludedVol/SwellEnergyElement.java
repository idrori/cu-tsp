package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.EnergyElement;
import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomList;

/**
 * Created by chen on 03/04/2016.
 */
public class SwellEnergyElement extends NonBondedEnergyElement {
    public SwellEnergyElement(DistanceMatrix distanceMatrix) {
        atoms = new AtomList(distanceMatrix.molecularSystem);
    }
    public void setAtoms() {};
    public double evaluate() {return 0;} ;
    public void set(Object o) {}

}
