package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Utils;
import meshi.util.filters.Filter;

/**
 * Created by chen on 03/04/2016.
 */
public class Swell extends ExcludedVol{
    public Swell(DistanceMatrix distanceMatrix,
                 EnergyInfoElement info) {
        super(distanceMatrix, null, 0, info, 0, null);
        comment = "Swell";
        energyElement = new SwellEnergyElement(distanceMatrix);
    }
    public EnergyInfoElement evaluate() throws EvaluationException {
        Atom atom1,atom2;
        double dEdD;
        double dDx,dDy,dDz;
        if (!on) {
            info.setValue(0);
            return info;
        }
        double energy = 0;
        DistanceLists nonBondedList = distanceMatrix.nonBondedList();
        for (DistanceList distanceRow : nonBondedList) {
            atom1 = distanceRow.atomOne.atom;
            for (Distance distance : distanceRow) {
                dDx = distance.dDistanceDx();
                dDy = distance.dDistanceDy();
                dDz = distance.dDistanceDz();
                atom2 = distance.atom2();
                if (distance.distance > distanceMatrix.rMax()) continue;
                double diff = distance.distance - distanceMatrix.rMax();
                energy += weight*diff * diff;
                dEdD = weight*2 * diff;
                atom1.addToFx(-dEdD * dDx);
                atom1.addToFy(-dEdD * dDy);
                atom1.addToFz(-dEdD * dDz);
                atom2.addToFx(dEdD * dDx);
                atom2.addToFy(dEdD * dDy);
                atom2.addToFz(dEdD * dDz);
            }
        }
        info.setValue(energy);
        return info;
    }//evaluate

}
