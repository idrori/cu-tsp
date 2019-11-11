package meshi.energy.hydrogenBond;

import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.util.UpdateableException;
import meshi.util.Utils;

import java.util.Iterator;

/**
 * @author amilev
 *
 * This class claculate the HydrogenBonds Energy over all the unbounded Hydrogen-Oxygen Pairs in the protein. 
 * To do so we keep an updatable-pairs-list with all the unbounded Hydrogen-Oxygen Pairs.
 */
public class HydrogenBondsEnergy extends NonBondedEnergyTerm {
    public static boolean debug = false;

    //--------------------------------- data fields ----------------------------------

    //*
    // * Accept a Distance between getAtoms if one of them is backbond Hydrogen and the other is backbond Oxygen
    // */
    //protected DistanceMatrix distanceMatrix;
    //protected Distance elementsList;
    protected HBondList hBondList;
    protected HydrogenBondsEnergyElement energyElement;
    private HydrogenBondsParametersList parametersList;
    private final double width;
   // private DistanceLists specialDisatnces = null;
    /**
     * @return Returns the distanceMatrix.
     */
    public final DistanceMatrix getDistanceMatrix() {
        return distanceMatrix;
    }

    public HBondList hBondList() {
        return hBondList;
    }
    /**
     * @return Returns the parametersList.
     */
    public final HydrogenBondsParametersList getParametersList() {
        return parametersList;
    }

    public int evalCounter = 0;
    public int elementEvalCounter =0;
    public int maxListSize = 0;
    public int maxFilterList = 0;
    
    //------------------------------------ constractor ---------------------------------------------

    public HydrogenBondsEnergy(){
        super();
        width = -1;
        energyElement = new HydrogenBondsEnergyElement(weight, width,distanceMatrix);
        throw new RuntimeException("Please do not use this constructor");
    }


    public HydrogenBondsEnergy (DistanceMatrix distanceMatrix,
                                HydrogenBondsParametersList parametersList,
                                EnergyInfoElement info,
                               HBondList hBondList, double width) {
        super(toArray(distanceMatrix,hBondList), //create the array of apdateableResources
               info,
              distanceMatrix);
        comment = "HydrogenBondsEnergy";
        this.distanceMatrix = distanceMatrix;
        this.hBondList = hBondList;
        energyElement = new HydrogenBondsEnergyElement(weight, width,distanceMatrix);
        this.parametersList = parametersList;
        this.width = width;

    }

   public HydrogenBondsEnergy (DistanceMatrix distanceMatrix,
                                HydrogenBondsParametersList parametersList,
                                EnergyInfoElement info,
                                HBondList hBondList,
                               DistanceLists specialDistances,
                               double width) {
        this(distanceMatrix, parametersList, info, hBondList, width);
    //   this.specialDisatnces = specialDistances ;

   }


    //-------------------------------------- methods ---------------------------------------
    public EnergyInfoElement evaluate() {
	if (! on) {
        info.setValue(0.0);
        return info;
    }
        evalCounter ++;
        double energy = 0, energyCurrent;
        if(hBondList.size() > maxListSize )
            maxListSize =   hBondList.size() ;

        int counter =0;
        for(Distance pair :hBondList) {
            counter++;
            energyElement.set(pair);
            energyCurrent = energyElement.evaluate();
            elementEvalCounter++;
            energy += energyCurrent;
            energyElement.freeElement();
        }
        if(counter > maxFilterList )
            maxFilterList = counter;

         if(! (energy >0 | energy <= 0) )
             throw new RuntimeException("Very weird");
         info.setValue(energy);
         return info;

    }

    public void evaluateAtoms() throws EvaluationException{
       if (! on) return;
       for (Distance pair : hBondList) {
            energyElement.set(pair);
            energyElement.evaluateAtoms();
            energyElement.freeElement();
        }
    }

    public void test(TotalEnergy energy, Atom atom) {
        debug = true;
        System.out.println("Testing " + this + " with " + atom + "\nWidth = " + width);
            energy.evaluate();
           DistanceLists nonBondedList = distanceMatrix.nonBondedList();
           for (DistanceList nonBondedRow : nonBondedList) {
               for (Distance nonBonded : nonBondedRow) {
                   if (nonBonded == null) {
                       System.out.println("Weird nonBondedList:");
                       for (Distance d : nonBondedRow) System.out.println(d);
                       throw new RuntimeException("Null distance in NonBondedList");
                   }
                   if ((nonBonded.atom1() == null) || (nonBonded.atom2() == null)) {
                       //System.out.println("Null distance in NonBondedList"+nonBonded);
                       System.out.println("Weird nonBondedList:");
                       for (Distance d : nonBondedRow) System.out.println(d);
                       throw new RuntimeException("Null atom in distance " + nonBonded);
                   }
                   if (nonBonded == null) throw new RuntimeException("nonBonded is null");

                   HB_DistanceAttribute hb_Attribute = (HB_DistanceAttribute) nonBonded.getAttribute(HB_DistanceAttribute.key);
                   if (hb_Attribute == null) continue;
                   HB_AtomAttribute atom1_attribute = (HB_AtomAttribute) nonBonded.atom1.atom.getAttribute(HB_AtomAttribute.key);
                   if (atom1_attribute == null) continue;
                   if ((atom != nonBonded.atom1()) && (atom != nonBonded.atom2())) continue;
                   System.out.println("Within " + this + " testing element " + nonBonded);
                   energyElement.set(nonBonded);
                   energyElement.test(energy, atom);
               }
           }
       }


}
