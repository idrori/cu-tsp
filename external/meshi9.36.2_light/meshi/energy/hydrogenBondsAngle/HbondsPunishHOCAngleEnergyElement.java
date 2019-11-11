package meshi.energy.hydrogenBondsAngle;

import meshi.geometry.*;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 14/03/2006
 * Time: 14:31:36
 * This class punish HydrogenBonds with HOC angles < 150
 * The energy function is zero when the angle is >=150 or when the distance between the Hydrogen
 * and the Oxygen is bigger then 3.5 A.
 */
public class HbondsPunishHOCAngleEnergyElement extends AbstractPunishAngleEnergyElement {
    //----------------------------------- constructors --------------------------------------

       public HbondsPunishHOCAngleEnergyElement(DistanceMatrix distanceMatrix,double weight)
       {
          super(distanceMatrix, weight);
       }

       public HbondsPunishHOCAngleEnergyElement(DistanceMatrix distanceMatrix,double weight,double xMax)
       {
          super(distanceMatrix, weight,xMax);
       }

    public HbondsPunishHOCAngleEnergyElement(DistanceMatrix distanceMatrix,double weight,double xMax, double maxAngle)
       {
          super(distanceMatrix, weight,xMax, maxAngle);
       }

       //---------------------------------- methods --------------------------------------------

       public String comment(){
           return "HbondsPunishHOCAngleEnergyElement";
       }

       public void setTheirdAtom(){
           theirdAtom = oAtom.residue().carbonylC() ;
       }

       public void setDAngleEnergyDatoms(double dAngleEnergyDAngle){
                    //dAngleEnergyDAngle = 2*a2*angleValue+b2;
                    dAngleEnergyDxHAtom = dAngleEnergyDAngle*angle.dangleDx1();//hAtom = atom1
                    dAngleEnergyDyHAtom = dAngleEnergyDAngle*angle.dangleDy1();
                    dAngleEnergyDzHAtom = dAngleEnergyDAngle*angle.dangleDz1();
                    dAngleEnergyDxOAtom = dAngleEnergyDAngle*angle.dangleDx2();//oAtom = atom2
                    dAngleEnergyDyOAtom = dAngleEnergyDAngle*angle.dangleDy2();
                    dAngleEnergyDzOAtom = dAngleEnergyDAngle*angle.dangleDz2();
                    dAngleEnergyDxTheidAtom = dAngleEnergyDAngle*angle.dangleDx3();//cAtom = atom3
                    dAngleEnergyDyTheidAtom = dAngleEnergyDAngle*angle.dangleDy3();
                    dAngleEnergyDzTheidAtom = dAngleEnergyDAngle*angle.dangleDz3();
       }

       public boolean setAngle() {
           Distance bond;
           if (hAtom.number() > oAtom.number()) bond = distanceMatrix.distance(hAtom,oAtom);
           else bond = distanceMatrix.distance(oAtom,hAtom);
           if (bond == null) return false;
           try {
                  angle = new QuickAndDirtyAngle(hAtom,oAtom,theirdAtom,distanceMatrix);
	       }
           catch(RuntimeException ex) {
	                   System.out.println("\n Failed to setResidue QuickAndDirtyAngle\n"+hAtom+"\n"+oAtom+"\n"+theirdAtom+"\n");
	                   System.out.println((new FreeDistance(hAtom,oAtom))+"\n"+(new FreeDistance(oAtom,theirdAtom)));
	                  throw ex;
	        }
            return true;
       }

    public void updateTheirdAtom(){
           updateTheirdAtom(weight);
    }
       public void updateTheirdAtom(double weight){
           if (! theirdAtom.frozen()) {
               theirdAtom.addToFx(-1 * deDxTheidAtom * weight); // force = -derivative
               theirdAtom.addToFy(-1 * deDyTheidAtom * weight); // force = -derivative
               theirdAtom.addToFz(-1 * deDzTheidAtom * weight); // force = -derivative
           }

       }

       public String toString() {
           return "HbondsPunishHOCAngleEnergyElement (h,o,c): "+hAtom.residueNumber()+
               " "+oAtom.residueNumber()+" "+theirdAtom.residueNumber()+" dis: "+distanceValue+"\n"+
               hAtom+"\n"+oAtom+"\n"+theirdAtom+"\nangle:"+angle+"\n"
               +"disEnergy "+distanceEnergy+" angleEnergy "+angleEnergy+" energy*weight "+energy*weight;
       }

}
