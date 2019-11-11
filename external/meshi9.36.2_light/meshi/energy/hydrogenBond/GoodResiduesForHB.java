package meshi.energy.hydrogenBond;


import meshi.geometry.Distance;
import meshi.geometry.DistanceLists;
import meshi.molecularElements.Residue;
import meshi.parameters.*;
import meshi.util.filters.Filter;
import meshi.molecularElements.atoms.*;

import java.util.Iterator;
public class GoodResiduesForHB implements Filter{
    public static final GoodSS goodSS = new GoodSS(); //Chen (22.9.11) made this object public static so that it cn be used outside
                                                      // this class.
    private IsHO isHO = new IsHO();
    private IsNot13 isNot13 = new IsNot13();
    private boolean firstTimeWornning = true;
    private DistanceLists specialDis = null;

    public   GoodResiduesForHB(){}
    public GoodResiduesForHB(DistanceLists disList){
        specialDis = disList;
    }

    public boolean accept(Object obj) {
        Distance distance = (Distance)obj;
        boolean ans =  isHO.accept(obj) && isNot13.accept(obj); // Chen (22.9.11 removed this condition && goodSS.accept(obj);
                                                                // This classification moves to the HB-level
        if (!ans) {
            return ans;
	}
        boolean frozen = distance.atom1() .frozen() & distance.atom2() .frozen() ;
         if (ans & frozen){
                if(firstTimeWornning )   //TODO check this - maybe frozen getAtoms should be looked at ?
                {
                    System.out.println("*****************************************************************************************************\n" +
                                                               "*** NOTE: frozen getAtoms does not consider as GoodResiduesForHB  !!! *****\n" +
                                                               "*****************************************************************************************************");
                    firstTimeWornning = false;
                }
                return false;
            }
        if(ans & specialDis != null && !specialDis.isEmpty() ){
            Distance pair;
            Iterator specDisIter = specialDis.iterator() ;
             Atom atom1,atom2;
            while( (pair  = (Distance) specDisIter .next()) != null){
                    atom1 = pair.atom1() ;
                    int atom1Num = atom1.number();
                    atom2 = pair.atom2() ;
                    int atom2Num = atom2.number();
                    if((distance.atom1.atom == atom1) | (distance.atom1.atom == atom2) | (distance.atom2.atom == atom1) | (distance.atom2.atom == atom2)){
	                        return (distance.atom1().number() == atom1Num & distance.atom2() .number() == atom2Num) |
                                                  (distance.atom1().number() == atom2Num & distance.atom2() .number() == atom1Num);
		    }
            }
        }
	return ans;

    }



    //--------------------------- internal class IsNot13 ----------------------------------

    /* this filter Accept Distance between Atoms just if there are at list 3 getAtoms away reset the sequense,
     *  (atomes that are less then 3 getAtoms away reset the sequense never create Hydrogen bonds with each other)
     */
    static class IsNot13 implements Filter{
        public boolean accept(Object obj) {
            Distance dis = (Distance)obj;
            if (!dis.atom1().chain().equalsIgnoreCase(dis.atom1().chain()))
                return true;
            int residueNumber1 = dis.atom1().residueNumber();
            int residueNumber2 = dis.atom2().residueNumber();
            return (!(Math.abs(residueNumber1-residueNumber2) < 3));
        }
    }


    //--------------------------- internal class GoodSS ---------------------------
            //chen 15.9.2011 preventing hydrogen bonds between residues of the same beta segment
    static class GoodSS implements Filter{
        public boolean accept(Object obj) {
            Distance dis = (Distance)obj;
            Atom atom1 = dis.atom1();
            Atom atom2 = dis.atom2();
            Residue residue1 = atom1.residue();
            Residue residue2 = atom2.residue();
            int number1 = residue1.number();
            int number2 = residue2.number();
            SecondaryStructure atom1SS = residue1.getSecondaryStructure();
            SecondaryStructure atom2SS = residue2.getSecondaryStructure();
            Residue high, low,current;
            if (atom1SS.equals(SecondaryStructure.COIL) || atom2SS.equals(SecondaryStructure.COIL))
                return false;
            if ((atom1SS.equals(SecondaryStructure.HELIX) || atom1SS .equals(SecondaryStructure.HELIX_OR_COIL)) &&
                (atom2SS.equals(SecondaryStructure.HELIX) || atom2SS .equals(SecondaryStructure.HELIX_OR_COIL) ))
                {
                    int residueNumber1 = dis.atom1().residueNumber();
                    int residueNumber2 = dis.atom2().residueNumber();
                    if (dis.atom1().type().backboneH() & dis.atom2().type().backboneO())
                         return ((residueNumber1-residueNumber2) == 4);
                    if (dis.atom2().type().backboneH() & dis.atom1().type().backboneO())
                         return ((residueNumber2-residueNumber1) == 4);
                    else
                       throw new RuntimeException("one of the getAtoms must be Hydrogen and the other Oxygen");

                } //todo check    helix
      //       return true;
            else {  // not a helix HB
                if ((atom1SS.equals(SecondaryStructure.SHEET) || atom1SS .equals(SecondaryStructure.SHEET_OR_COIL)) &&
                    (atom2SS.equals(SecondaryStructure.SHEET) || atom2SS .equals(SecondaryStructure.SHEET_OR_COIL))){
                    //We do not want it if both are reset the same beta-strand
                    if (Math.abs(number1-number2) < 20) { //longer beta strands are an exception
                        if (number1 < number2) {
                            low = residue1; high = residue2;
                        }
                        else low = residue2; high = residue1;
                        if (low.nextAtom() == null) return true; // It is assumed that missing residues ar in loops
                        current = low.nextAtom().residue();
                        while (current != high) {
                            if (current.getSecondaryStructure() != SecondaryStructure.SHEET) return true;
                            if (current.nextAtom() == null) return true;
                            current = current.nextAtom().residue();
                        }
                        return false; // they are reset the same strand
                    }
                    else return true; // probably not from the same strand ToDo find a better criterion
                }
                return false; //Not a beta HB either
            }

        }
    }
}
