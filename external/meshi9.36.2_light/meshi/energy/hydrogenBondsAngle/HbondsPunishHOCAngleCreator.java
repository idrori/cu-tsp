package meshi.energy.hydrogenBondsAngle;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 14/03/2006
 * Time: 14:47:17
 * To change this template use File | Settings | File Templates.
 */
public class HbondsPunishHOCAngleCreator extends EnergyCreator {

    HbondsPunishHOCAngleEnergy hBondsPunishHOCAngleEnergy;

    HydrogenBondsCreator hydrogenBondsCreator;


    /**
     * @param hydrogenBondsCreator
     **/
    public HbondsPunishHOCAngleCreator(HydrogenBondsCreator hydrogenBondsCreator) {
        super(InfoType.HB_PUNISH_HOC_ANGLE);
        this.hydrogenBondsCreator = hydrogenBondsCreator;
    }

    /**
     * @param weight
     * @param hydrogenBondsCreator
     **/

    /**
     * @param protein
     * @param distanceMatrix
     * @param commands - NOT IN USE !
     **/
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        HydrogenBondsEnergy HBE = hydrogenBondsCreator.getHydrogenBondsEnergy();
        double xMax;
        double maxAngle;
        try{
            xMax  = commands.firstWordFilter(HYDROGEN_BONDS_ANGLES ).secondWord(MAX_DISTANCE ).thirdWordDouble() ;
        }
        catch (RuntimeException r){
               xMax = -1;
        }
        try{
            maxAngle  = commands.firstWordFilter(HYDROGEN_BONDS_ANGLES ).secondWord(MAX_ANGLE ).thirdWordDouble() ;
        }
        catch (RuntimeException r){
               maxAngle = -1;
        }

        EnergyInfoElement info = new EnergyInfoElement(InfoType.HB_PUNISH_HOC_ANGLE, "HB angle HOC energy ", weight);

          if (xMax == -1){
                return term = new HbondsPunishHOCAngleEnergy(distanceMatrix,
                        HBE.hBondList(),
                        info,
                        hydrogenBondsCreator.getSpecialDis() ) ;
          }
          else {
              if (maxAngle == -1){
                return new HbondsPunishHOCAngleEnergy(distanceMatrix,
                          HBE.hBondList(),
                          info,
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax);
              }
              else
              return term = new HbondsPunishHOCAngleEnergy(distanceMatrix,
                          HBE.hBondList(),
                          info,
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax,
                        maxAngle);

          }
    }
}
