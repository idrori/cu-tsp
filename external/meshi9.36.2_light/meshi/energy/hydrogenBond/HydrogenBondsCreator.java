/*
 * Created on 17/11/2004
 * Window - Preferences - Java - Code Style - Code Templates
 */
package meshi.energy.hydrogenBond;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomCore;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

/**
 * @author amilev
 *
 * This class is used to create a HydroenBondEnergy.
 * The user does not need to know how to create such an energy term but only to give its weight. 
 */
public class HydrogenBondsCreator extends EnergyCreator implements KeyWords {

    //------------------------- data fields ---------------------------
    
    //private IsN isN = new IsN(); old version
    HydrogenBondsEnergy hydrogenBondsEnergy;
    protected HydrogenBondsParametersList parametersList= null;
    HBondList hBondList = null;

    /**
  * @return Returns the hydrogenBondsEnergy.
  */
 public final HydrogenBondsEnergy getHydrogenBondsEnergy() {	return hydrogenBondsEnergy;	}

    private DistanceLists specialDis = null;
    public final DistanceLists getSpecialDis(){return specialDis ;}

    int[ ] specialDisArray = null;
    public int[] getSpecialDisArray() {
        return specialDisArray;
    }

    private boolean antiParalel;
    public boolean getAntiParalel() {
        return antiParalel;
    }




    //------------------------- constructors ---------------------------

    public HydrogenBondsCreator() {
		super(InfoType.HYDROGEN_BONDS);
	}

    //25.9.11 Chen removed these apparently unused constructors and method
	/*public HydrogenBondsCreator(double weight) {
		super(weight);
	}

	public HydrogenBondsCreator(double weight , int[ ] specialDisArray,Protein protein,DistanceMatrix distanceMatrix,boolean antiParallel ){
		super(weight);
        readSpecialDistance(specialDisArray ,protein,distanceMatrix,antiParallel);
    }

     public HydrogenBondsCreator(int[ ] specialDisArray,Protein protein,DistanceMatrix distanceMatrix,boolean antiParallel){
		super(InfoType.HYDROGEN_BONDS);
        readSpecialDistance(specialDisArray ,protein,distanceMatrix,antiParallel);
    }

    private void readSpecialDistance(int[ ] specialDisArray,Protein protein,DistanceMatrix distanceMatrix, boolean antiParallel){
        this.antiParalel = antiParallel;
        if(specialDisArray .length % 2 != 0)
                        throw new MeshiException("HydrogenBondCreator Should get array of even length");
        this.specialDisArray = specialDisArray ;
        int length = specialDisArray.length;
             specialDis = new DistanceLists(100);
        if(antiParallel){
            for (int i = 0; i<length-1; i=i+2){
                specialDis.add(distanceMatrix.distance(protein.residue(specialDisArray [i]).getAtom("H"),protein.residue(specialDisArray [i+1]).getAtom("O") ) );
                specialDis.add(distanceMatrix.distance(protein.residue(specialDisArray [i+1]).getAtom("H"),protein.residue(specialDisArray [i]).getAtom("O") ) );
            }
        }
        else{
            for(int i=0;i<length-1;i=i+2){
                specialDis .add(distanceMatrix.distance(protein.residue(specialDisArray [i]).getAtom("O"), protein.residue(specialDisArray [i+1]+1).getAtom("H")));
                specialDis .add(distanceMatrix.distance(protein.residue(specialDisArray [i]).getAtom("H"), protein.residue(specialDisArray [i+1]- 1).getAtom("O") ) );
            }
        }
	Utils.print(specialDis);
    }
      */


    //------------------------------ methods ------------------------------

    /*
     * @param protein does not been used; must be given to implement an abstruct method of EnergyCreator.
     * @param distanceMatrix
     * @param commands gives the path to the data directory
     */
	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
                                           CommandList commands) {
        //create parameters
        if (parametersList== null)
            {
                //AtomList nitrogens = protein.getAtoms().filter(isN); old version
                parametersList = new HydrogenBondsParametersList(parametersDirectory(commands)+
                                                                 "/"+LENNARD_JONES_PARAMETERS);//,nitrogens);
            }
		for (AtomCore atom:distanceMatrix.molecularSystem){
		    if (atom.type().backboneH()) atom.atom.addAttribute(new HB_AtomAttribute(true,false));// changed 11/12/08
		    if (atom.type().backboneO()) atom.atom.addAttribute(new HB_AtomAttribute(false,true));// changed 11/12/08
        }
        double width = commands.firstWordFilter(HYDROGEN_BONDS).secondWord("width").thirdWordDouble();
                EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.HYDROGEN_BONDS, " Hydrogen bonds energy", weight);
                 
        if(specialDis != null ){

           hydrogenBondsEnergy = new HydrogenBondsEnergy(distanceMatrix,
                             parametersList,
							 energyInfoElement,
							 new HBondList(distanceMatrix, parametersList,specialDis ),
                    specialDis, width);
            if (hydrogenBondsEnergy.hBondList.size() == 0)
                throw new RuntimeException("This is weird 1");

        }
        else{
                hBondList =  new HBondList(distanceMatrix, parametersList);

            hydrogenBondsEnergy = new HydrogenBondsEnergy(distanceMatrix,
                                                      parametersList,
                                                      energyInfoElement,
                                                      hBondList,
                                                      width);
            if (hydrogenBondsEnergy.hBondList.size() == 0)
                throw new RuntimeException("This is weird 2, No hydrogen bonds.");
        }
            //Chen8.9.2011 bug fixed. The inputNewHBlist could be added a few times making some of the Hbonds counted many times including
        // the formation of weird "patterns" of a bond with itself.
        if (!distanceMatrix.energyTermsDistanceLists().contains(HBondList.inputNewHBList()))
                distanceMatrix.energyTermsDistanceLists() .add(HBondList.inputNewHBList() );
        return term = hydrogenBondsEnergy;
	}
	
	
    //---------------------------------------------------------------------------
	/* private static class IsN implements Filter {
       int[] sortBB_NITROGENS;
        
       //create a new SORT copy of BB_NITROGENS that are defined in interface AtomType
       public IsN(){             
       sortBB_NITROGENS = new int[BB_NITROGENS.length];
       for(int i = 0;i<BB_NITROGENS.length;i++)
       sortBB_NITROGENS[i] = BB_NITROGENS[i];
       Arrays.sort(sortBB_NITROGENS);
       }
        
       public boolean accept(Object obj) {
       return (Arrays.binarySearch(sortBB_NITROGENS,((Atom)obj).type) >= 0);
            
       }
       }//isN
    */
}
