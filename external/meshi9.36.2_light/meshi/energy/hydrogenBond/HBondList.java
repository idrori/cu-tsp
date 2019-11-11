package meshi.energy.hydrogenBond;

import meshi.geometry.*;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Updateable;
import meshi.util.filters.Filter;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * @author amilev
 *   Note: if the 2 getAtoms of the HB are frozen - we don't add them to the list.
 **/
public class HBondList extends ArrayList<Distance> implements Updateable {

    //--------------------------------- data fields -----------------------------
     public int insertionToList = 0;
    public int deletionFromList = 0;
    
    /*
     * The  PairsList updated every X steps.
     */
    private final int UPDATE_EVERY_X_STEPS = 50;
    ///*
    // * List of the candidate to hydrogen bonding given by the DistanceMatrix object
    // */
    //protected DistanceLists nonBondedList;
    /*
     * List of all the new HB elements that were added to the  hBondList in the last update call.
     */
    protected ArrayList<Distance> newhBondList;
    protected static HBdistanceLists inputNewHBList = new HBdistanceLists(new GoodResiduesForHB());
    /*
    * List of all the HB elements
    */
    //protected DistanceLists  hBondList;

    /*
     * holds the relevant row numbers
     */
    private ArrayList<Integer> relevantRows;
    /*
     * get the parameters (epsilon, sigma) of each element
     */
    private HydrogenBondsParametersList parametersList;
    private DistanceMatrix distanceMatrix;

    /*
     * Used to avoid reupdate in the same minimization step
     */
    private int numberOfUpdates =0;
    /*
     *  hBondList/this should be updated when countUpdates >= UPDATE_EVERY_X_STEPS
     */
    private int countUpdates = 50;
    /*
     * Filter for update
     */
    private GoodResiduesForHB goodResiduesForHB;
    private final IsAlive isAlive;

    /**
     * to find easily the last element that was in this list before the new update
     */
    private int oldSize = -1;

   public final int getOldSize() { return oldSize; }
    public final ArrayList<Distance> newhBondList() { return newhBondList;}
    public static HBdistanceLists inputNewHBList() { return inputNewHBList;}
    public final int countUpdates() {return countUpdates;}
    public final Filter filter;
    //-------------------------------- constructors --------------------------------

    /**
	 * @param distanceMatrix
	 */
	public HBondList(DistanceMatrix distanceMatrix,
                     HydrogenBondsParametersList parametersList) {
        super(100);
        goodResiduesForHB = new GoodResiduesForHB();
        this.parametersList = parametersList;
        this.distanceMatrix = distanceMatrix;
        //hBondList = new DistanceLists();
        newhBondList = new ArrayList<Distance>();
        relevantRows = new ArrayList<Integer>();
        //nonBondedList = distanceMatrix.nonBondedList();
        //update(nonBondedList);
        isAlive =  new IsAlive(distanceMatrix );
        update_dm();
        filter = new GoodResiduesForHB();
    }

    public HBondList(DistanceMatrix distanceMatrix,
                     HydrogenBondsParametersList parametersList,
                     DistanceLists specialDistances) {
        super();
        filter = new GoodResiduesForHB(specialDistances);
        goodResiduesForHB = new GoodResiduesForHB(specialDistances);
        this.parametersList = parametersList;
        this.distanceMatrix = distanceMatrix;
        relevantRows = new ArrayList<Integer>();
        isAlive =  new IsAlive(distanceMatrix );
        update_dm();
    }

    //--------------------------------------- methods ---------------------------------------

    public void print(){
            for (Distance pair:this) {
                System.out.println(pair);
            }

    }

    /* (non-Javadoc)
	 * @see meshi.util.Updateable#update(int)
	 */
	public void update(int numberOfUpdates) {
		if (numberOfUpdates == this.numberOfUpdates+1) {
		    this.numberOfUpdates++;
            //System.out.println("HBondLists: in update(int numberOfUpdates):");
		    update();
		}
		else if (numberOfUpdates != this.numberOfUpdates)
		    throw new RuntimeException("Something weird with HbondList.update(int numberOfUpdates)\n"+
                                       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
	}


    private void update() {
        oldSize = -1;
        newhBondList.clear();
        if (countUpdates == UPDATE_EVERY_X_STEPS) {//TODO can delete this if - i leave it for now for possible future needs
            int prevrousSize = size() ;
            updateTwoLists(inputNewHBList);
            int secondSize = size() ;
            countUpdates = 0;
            cleanList_dm();
            oldSize = size() - inputNewHBList .size() ;
            //System.out.println("HBondLists: prevouse size "+prevrousSize +" secondSize  "+secondSize +" current size "+size());
            if(size() > secondSize | secondSize < prevrousSize | size() < oldSize )
                throw new RuntimeException("HBondsList: problem at update - the list after adding elements is shorter or after cleaning is longer");
        }
        else if (countUpdates > UPDATE_EVERY_X_STEPS)
            throw new RuntimeException("Something weird with HbondList.update()\n"+
                                       "countUpdates = "+countUpdates+" UPDATE_EVERY_X_STEPS  = "+UPDATE_EVERY_X_STEPS);
        else {
            //System.out.println("HBondLists: in update: distanceMatrix.newHydrogenBondsList():");
            //distanceMatrix.newHydrogenBondsList().print();
            countUpdates++;
           updateTwoLists(inputNewHBList);
           cleanList_dm();
           oldSize = size() - inputNewHBList .size() ;
           if(size() < oldSize)
                     throw new RuntimeException("HBondLists: size < old size");
            }
           // test(1);
    }

	/**
	 * @param atomPairList of new distances that where added after  HbondList was updated in the last time
	 */
    private void updateTwoLists(HBdistanceLists atomPairList) {
        HB_DistanceAttribute distanceAttribute; // Chen 22.9.11 moved it here
        oldSize = -1;
        oldSize = size();
        for (DistanceList distanceList : atomPairList) {
            for (Distance pair:distanceList) {
                if (!isAlive.accept(pair))
                    System.out.println("HBondLists: there is a problem");
                //System.out.println("HBondLists: in updateTwoLists: pair: "+pair);
                //if (goodResiduesForHB.accept(pair)) {     this term is checked reset HBdistanceLists
                //System.out.println("!!!!!! HBondLists: in updateTwoLists: Goodpair !!!!! : "+pair);
                //Chen 22.9.11 Add all backbone hydrogen bonds to the list even those that we do not want to see in the
                // model
                if (GoodResiduesForHB.goodSS.accept(pair))
                    distanceAttribute = new HB_DistanceAttribute(true);
                else
                    distanceAttribute = new HB_DistanceAttribute(false);
                //
                distanceAttribute.setParameters((HydrogenBondsParameters) parametersList.parameters(pair));
                distanceAttribute.set(pair.atom1(), pair.atom2());
                pair.addAttribute(distanceAttribute);
                if ((pair.mode() == DistanceMode.INFINITE))
                    throw new RuntimeException("Weird pair in  updateTwoLists " + pair);
                add(pair);
                newhBondList.add(pair);
                insertionToList++;
                //}
            }
        }
    }

    /*
     * the update is done by first go over the heads getAtoms of the row in distance matrix
     * and just then (only if the head atom is Hydrogen or Oxygen) go over the rows themselves.
     * (more efficient ?)
     */
    private void update_dm(){
        // System.out.println("HBondLists: in update_dm: start");
        HB_DistanceAttribute distanceAttribute; //Chen 22.9.11 moved it here
        Atom headAtom;
        //System.out.println("HBondLists: dm  length = " +distanceMatrix .matrix.length );
        HB_AtomAttribute headAttribute;
        if (relevantRows.isEmpty()){
            // System.out.println("HBondList: relevantRows is empty");
            for (DistanceList nonBondeRow : distanceMatrix.nonBondedList()){
                //while((matrixRow = (MatrixRow) rows.next()) != null){
                    headAtom = nonBondeRow.atomOne.atom;
                    //System.out.println("HBondLists: headAtom = " +headAtom);
                    headAttribute = (HB_AtomAttribute) headAtom.getAttribute(HB_AtomAttribute.key);
                    if (headAttribute != null) {//means that it is H or O !
                        //System.out.println("HBondLists: headAttribute != null");
                        relevantRows.add(headAtom.number());
                        for (Distance pair : nonBondeRow) {
                            if (goodResiduesForHB.accept(pair)) {
                                if (!isAlive.accept(pair))
                                            throw new RuntimeException("need the second check !\n" +
                                                    pair + "\n" +
                                                    pair.atom1() + " " + pair.atom1().frozen() +
                                                    pair.atom2() + " " + pair.atom2().frozen());
                                //Chen 22.9.11 Add all backbone hydrogen bonds to the list even those that we do not want to see in the
                                // model
                                if (GoodResiduesForHB.goodSS.accept(pair)) {
                                            distanceAttribute = new HB_DistanceAttribute(true);
                                 } else {
                                            distanceAttribute = new HB_DistanceAttribute(false);
                                 }

                                 distanceAttribute.setParameters((HydrogenBondsParameters) parametersList.parameters(pair));
                                 distanceAttribute.set(pair.atom1(), pair.atom2());
                                 pair.addAttribute(distanceAttribute);
                                 if ((pair.mode() == DistanceMode.INFINITE))
                                      throw new RuntimeException("Weird pair in update_dm  " + pair);
                                 add(pair);
                                 insertionToList++;
                            }
                        }//while
                    }//if
            }
        }//if not empty
        else {
                for (Integer current:relevantRows){
                    DistanceList nonBondedRow = distanceMatrix.nonBondedList().search(current.intValue());
		            for (Distance pair:nonBondedRow) {
			            if(goodResiduesForHB.accept(pair)){
			                if(!isAlive .accept(pair ))
                                throw new RuntimeException("need the secod check !");
                            //Chen 22.9.11 Add all backbone hydrogen bonds to the list even those that we do not want to see in the
                            // model
				            if(GoodResiduesForHB.goodSS.accept(pair))
                                distanceAttribute = new HB_DistanceAttribute(true);
                            else
                                distanceAttribute = new HB_DistanceAttribute(false);
    			            distanceAttribute.setParameters((HydrogenBondsParameters) parametersList.parameters(pair));
                            //
			                pair.addAttribute(distanceAttribute);
			                // hBondList.add(pair);
			                if ((pair.mode() == DistanceMode.INFINITE)) throw new RuntimeException("Weird pair in update_dm #2  "+pair);
			                add(pair);
			                insertionToList ++;
			            }
		            }
                }//while
        }//else
        // System.out.println("HBondLists: in update_dm: end");
        // this.print();
    }



    /**
     * save only the live elements in this list
     */
    public void cleanList_dm(){
	Distance[] newInternalArray = new Distance[size()];
	int currentIndex =0;
	for(Distance pair:this){
	    if (isAlive .accept(pair )){
		newInternalArray [currentIndex] = pair ;
		currentIndex++;
               }
	    else
		deletionFromList++;
	}
	clear();
	for (Distance pair:newInternalArray) {
	    //if ((pair.mode() == DistanceMode.INFINITE)) throw new RuntimeException("cleanList_dm  "+pair);
	    if (pair != null) add(pair);
	}
    }

//       public Iterator withinRmaxIterator() {
//           return new WithinRmaxIterator(this);// WithinRmaxIterator is an internal class.
//       }


    //--------------------------- internal class IsWithInRmax ---------------------------

    static class IsWithInRmax implements Filter{
       private final double rMax;
        private DistanceMatrix dm;
        public IsWithInRmax(DistanceMatrix dm) {
            super();
            this.rMax = dm.rMax();
        }
        public boolean accept(Object obj) {
            Distance distance = (Distance)obj;
            //double dis = rMax - distance.distance();
            return ((rMax >= distance.distance()) & (!distance.dead()));
        }
    } //--------------------------- internal class IsAlive ---------------------------

    static class IsAlive implements Filter{
        DistanceMatrix dm;
        boolean firstTimeWornning = true;
        public IsAlive(DistanceMatrix  matrix){
            super();
            dm = matrix ;
        }

        public boolean accept(Object obj) {
            Distance distance = (Distance) obj;
            DistanceMode mode = distance.mode();
            if ((mode == DistanceMode.DEAD) || (mode == DistanceMode.INFINITE)) return false;
            return true;
        }
    }
    


    //--------------------------- internal class WithinRmaxIterator ---------------------------
    
    private  class WithinRmaxIterator implements Iterator   {
        IsWithInRmax isWithInRmax;
        DistanceLists list;
	int current;
	int listSize = 0;
        public WithinRmaxIterator(DistanceLists list){
        isWithInRmax = new IsWithInRmax(distanceMatrix);
	    this.list = list;
	    listSize = list.size();
	    current = 0;
        }

        /*
         * @return the next element that itas 2 getAtoms are within Rmax or Null if there is no such element
         */
        public Object next() {
	    if (current >= listSize) return null;
	    Object obj = list.get(current);
	    current++;
	    if (isWithInRmax.accept(obj)) return (obj);
	    else return next();
        }
	public boolean hasNext() {return (current < listSize);}
	public void remove() {throw new RuntimeException("do not do that");}
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }
}




