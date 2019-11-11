/*
 * Created on 26/01/2005
 * as part of meshi.1.5
 * 
 */
package meshi.energy.hydrogenBondsAngle;

import meshi.energy.EnergyInfoElement;
import meshi.energy.hydrogenBond.HBondList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
//import meshi.molecularElements.AtomPair;

/**
 * @author amilev
 *
 * This class punish HydrogenBonds with angles < 150 
 * The energy function is zero when the angle is >=150 or when the distance between the Hydrogen 
 * and the Oxygen is bigger then 3.5 A.
 */
public class HBondsPunishOHNAngleEnergy extends AbstractPunishAngleEnergy {
	
    public HBondsPunishOHNAngleEnergy(){}

    public HBondsPunishOHNAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      DistanceLists specialDis)
    {
        super(distanceMatrix,hBondList ,info,specialDis);
    }

    public HBondsPunishOHNAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      DistanceLists specialDis,
                                      double xMax)
    {
        super(distanceMatrix,hBondList ,info,specialDis,xMax);
    }

    public HBondsPunishOHNAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                       EnergyInfoElement info,
                                      DistanceLists specialDis,
                                      double xMax,
                                      double maxAngle)
    {
        super(distanceMatrix,hBondList ,info,specialDis,xMax,maxAngle);
    }


    /**
     * setResidue the energyElement to point reset the relavent instance
     */
    public void setEnergyElement() {
        energyElement = new HBondsPunishOHNAngleEnergyElement(distanceMatrix,weight);
    }

    public void setEnergyElement(double xMax) {
        energyElement = new HBondsPunishOHNAngleEnergyElement(distanceMatrix,weight,xMax);
    }

    public void setEnergyElement(double xMAx, double maxAngle) {
        energyElement = new HBondsPunishOHNAngleEnergyElement(distanceMatrix,weight,xMAx,maxAngle);
    }

    /**
     * setResidue the comment to be the relavent comment according to the instance
     */
    public void setComment() {
        comment = "HBondsPunishOHNAngleEnergy";
    }

}
