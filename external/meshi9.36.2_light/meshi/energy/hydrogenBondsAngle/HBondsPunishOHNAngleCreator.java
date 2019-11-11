/*
 * Created on 26/01/2005
 * as part of meshi.1.5
 * 
 */
package meshi.energy.hydrogenBondsAngle;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * @author amilev
  */
public class HBondsPunishOHNAngleCreator extends EnergyCreator {
	
	HBondsPunishOHNAngleEnergy hBondsPunishOHNAngleEnergy;
	HydrogenBondsCreator hydrogenBondsCreator;
	
	
	/**
	 * @param hydrogenBondsCreator
	 **/
	public HBondsPunishOHNAngleCreator(HydrogenBondsCreator hydrogenBondsCreator) {
		super(InfoType.HB_PUNIDH_OHN_ANGLE);
		this.hydrogenBondsCreator = hydrogenBondsCreator;
	}
	

	
	/**
	 * @param protein
	 * @param distanceMatrix
     * @param commands - NOT IN USE !
	 **/
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
		      CommandList commands)
	{
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

        EnergyInfoElement info = new EnergyInfoElement(InfoType.HB_PUNIDH_OHN_ANGLE, "HB angle HOC energy ", weight);
          if (xMax == -1){
        return new HBondsPunishOHNAngleEnergy(distanceMatrix,
                HBE.hBondList(),
                info,
                 hydrogenBondsCreator.getSpecialDis() ) ;
               }
          else {
              if (maxAngle == -1){
                   return new HBondsPunishOHNAngleEnergy(distanceMatrix,
                        HBE.hBondList(),
                        info,
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax) ;
              }
              else
              return new HBondsPunishOHNAngleEnergy(distanceMatrix,
                        HBE.hBondList(),
                        info,
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax,
                        maxAngle);

          }
    }

}
