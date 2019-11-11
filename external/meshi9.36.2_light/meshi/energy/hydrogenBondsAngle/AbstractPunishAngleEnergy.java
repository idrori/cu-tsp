package meshi.energy.hydrogenBondsAngle;

import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.hydrogenBond.HBondList;
import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;
import meshi.util.MeshiException;
import meshi.util.Utils;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 14/03/2006
 * Time: 14:13:56
 * To change this template use File | Settings | File Templates.
 */
public abstract class AbstractPunishAngleEnergy extends NonBondedEnergyTerm {
    protected HBondList hBondList;
    //private DistanceLists specialDisatnces = null;

    public AbstractPunishAngleEnergy(){}

    public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info)
    {
        super(toArray(distanceMatrix,hBondList),//updatabel resources
              info,
              distanceMatrix);
        this.hBondList = hBondList;
        setComment();
        setEnergyElement();
    }

    public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      double xMax)
    {
        super(toArray(distanceMatrix,hBondList),//updatabel resources
              info,
              distanceMatrix);
        this.hBondList = hBondList;
        setComment();
        setEnergyElement(xMax);
    }

     public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      double xMax,
                                      double maxAngle)
    {
        super(toArray(distanceMatrix,hBondList),//updatabel resources
              info,
              distanceMatrix);
        this.hBondList = hBondList;
        setComment();
        setEnergyElement(xMax,maxAngle);
    }

    public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      DistanceLists specialDis){
        this(distanceMatrix,hBondList,info);
        //specialDisatnces = specialDis;
    }

     public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      DistanceLists specialDis,
                                      double xMAx){
        this(distanceMatrix,hBondList,info,xMAx);
        //specialDisatnces = specialDis;
    }

     public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      EnergyInfoElement info,
                                      DistanceLists specialDis,
                                      double xMax,
                                      double maxAngle ){
        this(distanceMatrix,hBondList,info,xMax,maxAngle);
        //specialDisatnces = specialDis;
    }

    /**
     * setResidue the energyElement to point reset the relavent instance
     */
    public abstract void  setEnergyElement();

    /**
        * setResidue the energyElement to point reset the relavent instance
        */
       public abstract void  setEnergyElement(double xMax);

    /**
            * setResidue the energyElement to point reset the relavent instance
            */
           public abstract void  setEnergyElement(double xMax,double maxAngle);

    /**
     * setResidue the comment to be the relavent comment according to the instance
     */
    public abstract void  setComment();

    public EnergyInfoElement evaluate() throws EvaluationException {
        if (! on) {
            info.setValue(0.0);
            return info;
        }
        double energy = 0;
        double element_energy;


        for (Distance distance : hBondList) {
            energyElement.set(distance);

	        try {
		        element_energy = energyElement.evaluate();
	        }
	    catch (RuntimeException ex) {
		System.out.println("\n Failed to evaluate "+distance+"\n");
		throw ex;
	    }
            energy += element_energy;
            if (element_energy < 0) //energy should always be positive (if the angle is bad angle) or zero if it is good angle or pair of hidrogen-oxygen that are not connected
                throw new MeshiException(comment +": energy ("+energy +") should always be positive or zero: "+distance );

          //TODO add freeElement
        }
        info.setValue(energy);
        return info;
    }

    public void evaluateAtoms() throws EvaluationException{
        if (! on) return;
        for (Distance pair  : hBondList) {
            energyElement.set(pair);
            energyElement.evaluateAtoms();
            //TODO add freeElement
        }

    }

    public void test(TotalEnergy totalEnergy, Atom atom){
        System.out.println("Start Test "+comment);
        if (! on) System.out.println(""+this +" is off");
        for (Distance pair  : hBondList) {
            energyElement.set(pair);
            energyElement.test(totalEnergy,atom);
            //TODO add freeElement
        }
        System.out.println("End Test "+comment);
    }
}
