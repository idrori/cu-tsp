/*
 * Created on 22/8/2005
 *
 */
package meshi.energy.hydrogenBond;

import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.*;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.util.Utils;

/**
 * @author amilev
 *
 * This class is an energy element for an Hydrogen Bond. Two getAtoms are involved: Hydrogen atom (hAtom,
 * and Oxygen Atom (oAtom).
 * The Hydrogen Bond forces energy between hAtom and oAtom is similar to
 * Lennard-Jones in the range [minimum:infinity]
 * where the x value of the minimum is sqrSixOf2*Sigma.
 * in the range of [0:minimum] the energy function is decreasing very slowly.
 * E(Distance) =
 *   0.0001*Epsilon*(Distance - sqrSix(2)*Sigma)*(Distance - sqrSix(2)*Sigma) - Epsilon  //Distance < sqrSix(2)*Sigma
 *         [4*sigma^6*Epsilon*(sigma^6/Distance^12 - 1/Distance^6)]*[(rMax-Distance)^2/(rMax-Distance)^2+width
 *   (a graph shown this function can be found in file: UnderZeroLJEnergyElement.gif)
 */
public class HydrogenBondsEnergyElement extends NonBondedEnergyElement {

    //---------------------------------- data filds -----------------------------------

    /*
     * weight given to this energy in the TotalEnergy element
     */
    private double weight;  
    private double energy;
    protected Atom oAtom, hAtom;
    private double deDxOAtom, deDyOAtom, deDzOAtom;
    private double deDxHAtom, deDyHAtom, deDzHAtom;
    private double dEdD, dEdX, dEdY, dEdZ;
    
    private int hFactor,oFactor;
    private double HBfactor; //Chen 22.9.11 this factor is 1 for Hbonds within alpha-helix and between beta strands
                                 // Negative value for other less wanted bonds
    private double HBwidth; //Chen 22.9.11 user defined  for Hbonds within alpha-helix and between beta strands
                                                              // small value for other less wanted bonds
    public final double BAD_HB_FACTOR = -2;//chen 22.9.11
    public final double BAD_HB_WIDTH = 2;
    
    //--------- fields needed for energy calculationm -------

    public final double width;
    public static final double EPSILON = 0.01;
    public static final double MIN = 1.9;

    private double dis, disMinusMin, disMinusMax,   disMinusMin3, disMinusMax3,   disMinusMin4, disMinusMax4;
    private double disMinusMin4PlusWidth, disMinusMax4PlusEpsilon, disMinusMin4PlusWidth2, disMinusMax4PlusEpsilon2;
    private double rMax;
    private Atom atom1, atom2;
    /*
     * contact is used to make sure the function get to 0 in rMax
     */
    private double e1,e2,dE1dDis,dE2dDis;

    //true if no setResidue was done or explisit call was done to freeElement
    private boolean free = true;
    /*
       /*
     * An Element is built according to distance between two getAtoms
     */
    private Distance distance;
    /*
     * Attribute of distance contains:
     * is this relevant hydrogenBond;
     * its parameters;
     * its energy and derivative
     */
    public HB_DistanceAttribute hb_Attribute;
    //private LennardJonesParameters parameters;

    
  /*  Chen 23.9.11 Removed this constructor, which is not used any more
       public HydrogenBondsEnergyElement() {
          rMax = DistanceMatrix.rMax();
        width = -1;
        throw new RuntimeException("Please do not use this constructor");
      }  */

    public HydrogenBondsEnergyElement(double weight, double width,DistanceMatrix distanceMatrix) {
        this.weight = weight;
        rMax = distanceMatrix.rMax();
        this.width = width;
    }

    /*
    * @parm obj should be Distance
    */

    public void set(Object obj){
        double separation;

        if (width < 0) throw new RuntimeException("width = " + width);
        free = false;
        distance = (Distance) obj;
        atom1 = distance.atom1();
        atom2 = distance.atom2();
        atoms = new AtomList(AtomPair.ATOM_PAIR_CAPACITY,atom1.molecularSystem);
	    atoms.add(atom1);
	    atoms.add(atom2);
        HB_AtomAttribute atom1_attribute = (HB_AtomAttribute) atom1.getAttribute(HB_AtomAttribute.key);
        if (atom1_attribute.isH) {
            hFactor = oFactor = -1;    
            hAtom = atom1;
            oAtom = atom2;
        }
        else {
           hFactor = oFactor = 1;
            hAtom = atom2;
            oAtom = atom1;
        }
        //nAtom = hAtom.residue().amideN();    

        hb_Attribute = (HB_DistanceAttribute) distance.getAttribute(HB_DistanceAttribute.key);
        //update the parameters of the element:
        //parameters = hb_Attribute.parameters; 
        //parameters = (LennardJonesParameters) parametersList.parameters(distance);
        if (hb_Attribute.isHbond) {
            HBfactor = 1;
            HBwidth = width;
        }
        else {
            HBwidth = BAD_HB_WIDTH;
            separation = atom1.residueNumber()-atom2.residueNumber();
            if (Math.abs(separation) >= 6)
                HBfactor = BAD_HB_FACTOR;
            else HBfactor = 0.1;

        }
        dis                   = distance.distance();
        disMinusMin   = dis- MIN;
        disMinusMax = dis-rMax;
       disMinusMin3   = disMinusMin*disMinusMin*disMinusMin;
       disMinusMax3=  disMinusMax*disMinusMax*disMinusMax;
       disMinusMin4  =   disMinusMin3* disMinusMin;
       disMinusMax4=   disMinusMax3* disMinusMax;
        disMinusMin4PlusWidth = disMinusMin4 + HBwidth;
       disMinusMax4PlusEpsilon = disMinusMax4+EPSILON ;
        disMinusMin4PlusWidth2 = disMinusMin4PlusWidth * disMinusMin4PlusWidth;
      disMinusMax4PlusEpsilon2 = disMinusMax4PlusEpsilon*disMinusMax4PlusEpsilon;
//        if (hb_Attribute.isHbond) HBfactor = 1;
//        else HBfactor = BAD_HB_FACTOR;
        
    }    

       //22.9.11 Chen removed as apparently this method has no usage.
   /* public double evaluate(double weight) {
        energy = updateEnergy();
        updateAtoms(weight);
        return energy * weight;
    }*/


    public double evaluate(){
        energy = updateEnergy();
        updateAtoms();
//        if (energy * weight*HBfactor > 1)
//            throw new RuntimeException("\n"+atom1+"\n"+atom2+"\n"+dis+"\n"+energy+"  "+HBfactor);
//        if (HBfactor< 0)
//            System.out.println(distance+"\n"+energy);
        return energy * weight*HBfactor;
     }


    /*
    * free all the global vars to avoid missuse
    */
    public void freeElement(){
        free = true;
        distance = null;
        hb_Attribute = null;
        atoms = null;
        oAtom = hAtom = atom1 =atom2 =null;
        dEdD =dEdX =deDxHAtom =deDxOAtom =dEdY = 0;
        deDyOAtom =deDyHAtom =dEdZ = deDzHAtom =deDzOAtom =0;
        energy =0;
        hFactor =oFactor = 0;
    }
    
   /**
    * energy and dirivarives calculation.
    **/
   public double updateEnergy() {
       double dis = distance.distance();
       double rMaxMinusDis = rMax - dis;
       double energy1, dE1dD;
        e1 = disMinusMin4 / disMinusMin4PlusWidth - 1;
        dE1dDis = 4 * disMinusMin3 * (1 / disMinusMin4PlusWidth - disMinusMin4 / disMinusMin4PlusWidth2);
       e2           =  disMinusMax4/disMinusMax4PlusEpsilon;
       dE2dDis = 4*disMinusMax3*(1/disMinusMax4PlusEpsilon-disMinusMax4/disMinusMax4PlusEpsilon2);
       energy  = e1*e2;
       dEdD    = dE1dDis*e2+e1*dE2dDis;
        if (energy > 0)
            throw new RuntimeException("In hydrogenBondsEnergyElement: Energy should not get this values ! "+
                    atom1+"\n"+
                    atom2+"\n"+
                    dis+"\n" +
                    disMinusMin4 + " " + disMinusMin4PlusWidth + " " + e1 + "\n" +
                    e2+" "+energy);
        dEdX = dEdD*distance.dDistanceDx();
        dEdY = dEdD*distance.dDistanceDy();
        dEdZ = dEdD*distance.dDistanceDz();

        deDxOAtom =    dEdX * oFactor;
        deDyOAtom =    dEdY * oFactor;
        deDzOAtom =    dEdZ * oFactor;
        deDxHAtom = -1*dEdX * hFactor;
        deDyHAtom = -1*dEdY * hFactor;
        deDzHAtom = -1*dEdZ * hFactor;

// following fields of the variable pair is need for the unit hydrogenBondsPairs
        hb_Attribute.set(deDxOAtom,
			      deDyOAtom,
			      deDzOAtom,
			      deDxHAtom,
			      deDyHAtom,
			      deDzHAtom,
			      energy);
         return energy;
   }

    private void updateAtoms(){
        updateAtoms(weight*HBfactor);
    }

    private void updateAtoms(double weight){
       if (! oAtom.frozen()) {
           oAtom.addToFx(-1 * deDxOAtom * weight); // force = -derivative
           oAtom.addToFy(-1 * deDyOAtom * weight); // force = -derivative
           oAtom.addToFz(-1 * deDzOAtom * weight); // force = -derivative   
       }
       if (! hAtom.frozen()) {
           hAtom.addToFx(-1 * deDxHAtom * weight);
           hAtom.addToFy(-1 * deDyHAtom * weight);
           hAtom.addToFz(-1 * deDzHAtom * weight);
       }
// following fields of the variable pair is need for the unit hydrogenBondsPairs
        hb_Attribute.set (oAtom,
			      hAtom);
   }
	
   public final double deDxOAtom(){return deDxOAtom;}
   public final double deDyOAtom(){return deDyOAtom;}
   public final double deDzOAtom(){return deDzOAtom;}
   public final double deDxHAtom(){return deDxHAtom;}
   public final double deDyHAtom(){return deDyHAtom;}
   public final double deDzHAtom(){return deDzHAtom;}
    //    public UnderZeroLJEnergyElement ljElement(){return LJElement;} //  
   public final double energy(){return energy;}
   public final Atom oAtom() {return oAtom;}
   public final Atom hAtom() {return hAtom;}  
             
   public String toString() {
       if ((oAtom == null) & (hAtom == null)) return "HydrogenBondEnergyElement: - getAtoms not setResidue yet";
       if ((oAtom == null) | (hAtom == null)) throw new RuntimeException("HydrogenBondEnergyElement: This is weird\n"+
                                                                         "hAtom = "+hAtom+"\n"+
                                                                         "oAtom = "+oAtom); 
       if (free) return "HydrogenBondsEnergyElement was free (no values)!";
       return "HydrogenBondsEnergyElement (h,o): "+hAtom.residueNumber()+
           " "+oAtom.residueNumber()+" dis: "+distanceValue()+"\n"+
						 hAtom+"\n"+oAtom;
       //return "HydrogenBondsEnergyElement";
   }

    public final double distanceValue() {
        return distance.distance();
    }
    protected void setAtoms(){
	throw new RuntimeException("setAtoms() may not be used by HydrogenBondsEnergyElement for "+
				   "efficiency.");
    }
    
}
