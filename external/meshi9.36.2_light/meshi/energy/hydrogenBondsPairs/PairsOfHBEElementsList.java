package meshi.energy.hydrogenBondsPairs;

import meshi.energy.hydrogenBond.HB_DistanceAttribute;
import meshi.energy.hydrogenBond.HBondList;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.geometry.Distance;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.util.Updateable;
import meshi.util.filters.Filter;

import java.util.ArrayList;
import java.util.Iterator;


/*
 * All the elements in this list are pairs of HB that are good pairs
 */
public class PairsOfHBEElementsList extends ArrayList<PairOfHydrogenBondsElements> implements Updateable {

    //--------------------------------------- data -------------------------------------------

    public int insertion =0;
    public int deletions =0;
     /*
     * update list of the hydrogen bonds, used to update this list whenever the hBondList has been changed  
     */
    protected HBondList hBondList;
    private int numberOfUpdates;
    private boolean debug = false;
    private DistanceMatrix distanceMatrix;

    public final HBondList hBondList(){return hBondList;}

    private IsWithInRmax isWithinRmax;
    //-------------------------------------- constructors ----------------------------------------
    
    public PairsOfHBEElementsList(){
        super();
    }

    public PairsOfHBEElementsList(HydrogenBondsEnergy hbe, DistanceMatrix distanceMatrix){
        super();
        hBondList = hbe.hBondList();
        createPairsList(hBondList );
        isWithinRmax = new IsWithInRmax(distanceMatrix);
        this.distanceMatrix = distanceMatrix;
    }

     //-------------------------------------- methods ----------------------------------------------
    
    public void print()
    {
        for (PairOfHydrogenBondsElements pair:this)
        {
            System.out.println(pair);
        }

    }

    public void print(Filter filter )
    {
        for (PairOfHydrogenBondsElements pair:this)
        {
            if(filter.accept(pair))
                System.out.println(pair);
        }

    }




    public void update(int numberOfUpdates) {
        if (numberOfUpdates == this.numberOfUpdates+1) {
            update();
            this.numberOfUpdates++;
        }
        else if (numberOfUpdates != this.numberOfUpdates) 
            throw new RuntimeException("Something weird with PairsOfHBEElementsList.update(int numberOfUpdates)\n"+
                                       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }
	    
    public void update() {
        if (hBondList.countUpdates() == 0){ //when countUpdates == 0 it meens that hBondList was restart
            updateNew(hBondList.newhBondList(),hBondList .getOldSize() );
            cleanList();
            if(debug){
                System.out.println("HBE:update:in first if:el: "+ size());
            }
        }
      else
            updateNew(hBondList.newhBondList(),hBondList .getOldSize() );
      cleanList();
    }

    private void createPairsList(ArrayList<Distance> hBonds){
        HB_DistanceAttribute attribute; // Chen 23.9.11
    	int length = hBonds.size();

        if (hBonds == null)
            throw new RuntimeException("hBonds is null");
        if (hBonds.get(0) == null)
            throw new RuntimeException("hBonds.get(0) is null\n hBonds = "+hBonds+"\n size is "+hBonds.size());
        if (hBonds.get(0).atom1() == null)
            throw new RuntimeException("hBonds.get(0).atom1() is null");

        PairOfHydrogenBondsElements setPair = new PairOfHydrogenBondsElements(hBonds.get(0).atom1().molecularSystem);
        Distance pair1 , pair2;
    	for(int i=0;i<length;i++){
            pair1 = hBonds.get(i);
            attribute = (HB_DistanceAttribute) pair1.getAttribute(HB_DistanceAttribute.key); //Chen 23.9.11
            if (attribute.isHbond) {                                                         //Chen 23.9.11
    		    for(int j=i+1;j<length;j++){
                        pair2 = hBonds.get(j);
                        attribute = (HB_DistanceAttribute) pair2.getAttribute(HB_DistanceAttribute.key); //Chen 23.9.11
                        if (attribute.isHbond) {                                                         //Chen 23.9.11
                            setPair.set(pair1,pair2);
                            if (setPair.lookAtThisPair) {
                                PairOfHydrogenBondsElements pairOfHydrogenBondsElements = new PairOfHydrogenBondsElements(setPair);
                                add(pairOfHydrogenBondsElements);
                                insertion ++;
                            }
                        }
                    }
                }
        }
    }
  
    public void updateNew(ArrayList<Distance> newlist,int lastIndexBeforeNewElements){
        if (lastIndexBeforeNewElements == -1)
                throw new RuntimeException("problem ... this method should be called just if HBondLists called to updateTwoLists !");
        if (newlist .size() !=0) {
            PairOfHydrogenBondsElements setPair = new PairOfHydrogenBondsElements(newlist.get(0).atom1().molecularSystem);
            Distance pair1 , pair2;
            for(int j=0;j<newlist.size();j++){
                for(int i=0;i<lastIndexBeforeNewElements ;i++){
                     pair1 = hBondList .get(i);
                    pair2 = newlist.get(j);
                    setPair.set(pair1, pair2);
                    if(setPair.lookAtThisPair) {
                        PairOfHydrogenBondsElements pairOfHydrogenBondsElements = new PairOfHydrogenBondsElements(setPair);
                        add(pairOfHydrogenBondsElements);
                        insertion++;
                    }
                }//while
            }//for
            createPairsList(newlist);
        }//if
    }

     /**
        * save only the live elements inthis list
        */
       public void cleanList(){
           Object[] newInternalArray = new Object[size() ];
           int currentIndex =0;
           for (PairOfHydrogenBondsElements pair:this){
	       if ((pair != null) && (isWithinRmax .accept(pair ))){
                    newInternalArray [currentIndex] = pair ;
                   currentIndex++;
                }
               else
                    deletions ++;
           }
           clear();
	   for (Object obj:newInternalArray)
	       add((PairOfHydrogenBondsElements)obj);
       }

    public Iterator withinRmaxPairsIterator() {
        return new WithinRmaxPairsIterator(this,distanceMatrix);
    } 

    //---------------------------------------------------------------------------------------
    

    //--------------------------- internal class IsWithInRmax ---------------------------

    static class IsWithInRmax implements Filter{
        private double rMax;
        public IsWithInRmax(DistanceMatrix dm){
            super();
            rMax = dm.rMax() ;
        }

        public boolean accept(Object obj) {
	    if (obj == null) return false;
            PairOfHydrogenBondsElements pair = (PairOfHydrogenBondsElements)obj;
            return !(pair.HOelement1.dead() || pair.HOelement2.dead()) && (rMax >= pair.HOelement1.distance()) &
                    (rMax >= pair.HOelement2.distance());
        }
    }
    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }
    //---------------------------------------------------------------------------------------
    
    //----------------------------- parivate class WithinRmaxPairsIterator ----------- ------------

    //---------------------------------------------------------------------------------------
    private  class WithinRmaxPairsIterator  implements Iterator{

        private IsWithInRmax isWithInRmax;
	private int current;
	private ArrayList list;
	private int listSize;
        public WithinRmaxPairsIterator(ArrayList list,DistanceMatrix distanceMatrix) {
            super();
	    this.list = list;
	    listSize = list.size();
            isWithInRmax = new IsWithInRmax(distanceMatrix);
	    current = 0;
        }
        
        /*
         * @return the next element that in each of its 2 pairs the 2 getAtoms are within Rmax or Null if there is no such element
         */
        public Object next() {
	    while (current < listSize) {
		Object currnetObject = list.get(current);
		current++;
		if (isWithInRmax.accept(currnetObject)) return currnetObject;
	    }
	    return null;
	}
	public boolean hasNext() {return (current<listSize);}
	public void remove() {throw new RuntimeException("do not do that");}
    }
}
