package meshi.energy.hydrogenBondsAngle;

import meshi.energy.EnergyInfoElement;
import meshi.energy.hydrogenBond.HBondList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 14/03/2006
 * Time: 14:43:32
 * This class punish HydrogenBonds with HOC angles < 150 
 * The energy function is zero when the angle is >=150 or when the distance between the Hydrogen
 * and the Oxygen is bigger then 3.5 A.
 */
public class HbondsPunishHOCAngleEnergy extends AbstractPunishAngleEnergy {

    public HbondsPunishHOCAngleEnergy(){}

    public HbondsPunishHOCAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      DistanceLists specialDis)
    {
        super(distanceMatrix,hBondList ,info,specialDis);
    }

     public HbondsPunishHOCAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      DistanceLists specialDis,
                                      double xMax)
    {
        super(distanceMatrix,hBondList ,info,specialDis,xMax);
    }

      public HbondsPunishHOCAngleEnergy(DistanceMatrix distanceMatrix,
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
        energyElement = new HbondsPunishHOCAngleEnergyElement(distanceMatrix,weight);
    }

    /**
     * setResidue the energyElement to point reset the relavent instance
     */
    public void setEnergyElement(double xMax) {
        energyElement = new HbondsPunishHOCAngleEnergyElement(distanceMatrix,weight,xMax);
    }


      /**
     * setResidue the energyElement to point reset the relavent instance
     */
    public void setEnergyElement(double xMax,double maxAngle) {
        energyElement = new HbondsPunishHOCAngleEnergyElement(distanceMatrix,weight,xMax,maxAngle);
    }
    /**
     * setResidue the comment to be the relavent comment according to the instance
     */
    public void setComment() {
        comment = "HbondsPunishHOCAngleEnergy";
    }

}
