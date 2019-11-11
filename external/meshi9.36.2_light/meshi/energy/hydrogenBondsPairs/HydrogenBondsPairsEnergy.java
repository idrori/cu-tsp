
/*
 * Created reset 26/12/2004
 *
 */
package meshi.energy.hydrogenBondsPairs;

import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.hydrogenBond.HB_DistanceAttribute;
import meshi.energy.hydrogenBond.SynchronizedUpdateException;
import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;

import java.util.Iterator;

/**
 * @author amilev
 */
public class HydrogenBondsPairsEnergy extends NonBondedEnergyTerm {
    private final double width;

    //--------------------------------- data ------------------------------------------

    public PairsOfHBEElementsList getPairsOfHBEElementsList() {
        return pairsOfHBEElementsList;
    }

    //protected HydrogenBondsPairsEnergyElement energyElement;
    protected PairsOfHBEElementsList pairsOfHBEElementsList;


    private static final double DEFAULT_PUNISHMENT = 10;        //TODO check this parameter
    private static final double DEFAULT_HPUNISHMENT = 30;    //TODO check this parameter

    public int evalCounter = 0;
    public int elementEvalCounter = 0;
    public int maxListSize = 0;
    public int maxFilterListSize = 0;
    public HelixParametersList helixParametersList = null;
    public BetaParametersList betaParametersList = null;
    public HelixParametersList helixBackgroundParametersList = null;
    public BetaParametersList betaBackgroundParametersList = null;
    //--------------------------------- constructors ----------------------------------
    
    /**
     * 
     */
    public HydrogenBondsPairsEnergy() {
        super();	
        width = -1;
    }


    public HydrogenBondsPairsEnergy(DistanceMatrix distanceMatrix,
                                    HelixParametersList helixParametersList,
                                    BetaParametersList betaParametersList,
                                    HelixParametersList helixBackgroundParametersList,
                                    BetaParametersList betaBackgroundParametersList,
                                    PairsOfHBEElementsList pairsOfHBEElementsList,
                                     EnergyInfoElement info,
                                    int [] specialDisArray,
                                    boolean antiParallel, double width) {
        this(distanceMatrix,
             helixParametersList,
             betaParametersList,
             helixBackgroundParametersList,
             betaBackgroundParametersList,
             pairsOfHBEElementsList,
             info,
             DEFAULT_PUNISHMENT,
             DEFAULT_HPUNISHMENT,
             specialDisArray,antiParallel,width);
    }

    /**
     * @param distanceMatrix  the distance matrix of the protein
     * @param helixParametersList parameters for pairs of two HB in helices
     * @param betaParametersList   parameters for pairs of two HB in betas
     * @param helixBackgroundParametersList   parameters for pairs of candidate HB in helicesin which at least one candidate is exsualy HB
     * @param betaBackgroundParametersList   parameters for pairs of candidate HB in betasin which at least one candidate is exsualy HB
     * @param pairsOfHBEElementsList   list of all pairs that the term is worked reset
     * @param info  of the term
     * @param punish  for unseen patterns
     * @param hpunish for patterns of HB where the H getAtoms is the same (two bounds reset the same H atom)
     * @param specialDisArray array of pairs of residues that are known to have HB in a beta sheet
     * @param antiParallel true if the residues in  specialDisArray are in parallel configuration or flase if in anti-parallel
     */
   public HydrogenBondsPairsEnergy(DistanceMatrix distanceMatrix,
                                    HelixParametersList helixParametersList,
                                    BetaParametersList betaParametersList,
                                    HelixParametersList helixBackgroundParametersList,
                                    BetaParametersList betaBackgroundParametersList,
                                    PairsOfHBEElementsList pairsOfHBEElementsList, 
                                    EnergyInfoElement info,
                                    double punish,
                                    double hpunish,
                                    int[] specialDisArray,
                                    boolean antiParallel,
                                    double width){
        
        super(toArray(distanceMatrix,pairsOfHBEElementsList.hBondList(),pairsOfHBEElementsList),
              info,
              distanceMatrix);
        comment = "HB-PairsEnergy";
        this.pairsOfHBEElementsList = pairsOfHBEElementsList;
        energyElement = new HydrogenBondsPairsEnergyElement(helixParametersList,
                                                            betaParametersList,
                                                            helixBackgroundParametersList,
                                                            betaBackgroundParametersList,
                                                            weight,punish,hpunish,
                                                            specialDisArray,antiParallel,width,distanceMatrix);
        this.helixParametersList= helixParametersList;
        this.betaParametersList=betaParametersList;
        this.helixBackgroundParametersList= helixBackgroundParametersList;
        this.betaBackgroundParametersList=betaBackgroundParametersList;
        this.width = width;
    }
	

    //---------------------------------------- methods ------------------------------------------
    
    /* (non-Javadoc)
     * @see meshi.energy.AbstractEnergy#evaluate()
     */
    public EnergyInfoElement evaluate()throws EvaluationException {
        if (! on) {
            info.setValue(0.0);
            return info;
        }

        double energy = 0;
        Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();
        PairOfHydrogenBondsElements pair;

        while ((pair  = (PairOfHydrogenBondsElements) pairsIter.next()) != null) {

            if (pair.isHelixPair() | pair.isBetaPair()) //TODO change when the system can deal with "ALL" types
            {
                ((HydrogenBondsPairsEnergyElement)energyElement).set(pair);
                if ( ((HydrogenBondsPairsEnergyElement)energyElement).factor > 10000)
                    System.out.println(((HydrogenBondsPairsEnergyElement)energyElement).factor);
                energy += energyElement.evaluate();
                ((HydrogenBondsPairsEnergyElement )energyElement) .freeElenet();
            }
        }
        info.setValue(energy);
        return info;
    }

    /* (non-Javadoc)
     * @see meshi.energy.AbstractEnergy#evaluateAtoms()
     */
    public void evaluateAtoms() throws EvaluationException{
        if (! on) return;
        Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();
        PairOfHydrogenBondsElements pair;
        while ((pair  = (PairOfHydrogenBondsElements) pairsIter.next()) != null) {
            {
                 if (pair.isHelixPair() | pair.isBetaPair()){ //TODO change when the system can deal with "ALL" types
                       energyElement.set(pair);
                       energyElement.evaluateAtoms();
                       ((HydrogenBondsPairsEnergyElement )energyElement) .freeElenet();
                 }
            }
        }
    }

    /* (non-Javadoc)
     * @see meshi.energy.AbstractEnergy#test(meshi.energy.TotalEnergy, boolean, boolean)
     */
    public void test(TotalEnergy totalEnergy,Atom atom){
        System.out.println("Start Test "+comment);
        System.out.println("atom = "+atom);
         System.out.println("hBondList: ");
        for(int i =0;i<pairsOfHBEElementsList .hBondList .size();i++){
                    Distance currentDis = pairsOfHBEElementsList.hBondList .get(i);
                    if ((currentDis.atom1.atom == atom) | (currentDis.atom2.atom == atom) )
                        System.out.println(currentDis);
        }

        System.out.println("pairsOfHBEElementsList: ");
        PairOfHydrogenBondsElements pair;
        for(Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();pairsIter.hasNext();){
            pair = (PairOfHydrogenBondsElements)pairsIter.next();
            if (pair == null) {
                  System.out.println("Null pair");
                   continue;
            }
            try {
            if((pair.HOelement1.atom1.atom == atom) | (pair.HOelement2 .atom1.atom == atom) | (pair.HOelement1.atom2.atom == atom) |
	       (pair.HOelement2.atom2.atom == atom) )
                        System.out.println("pair is "+pair);
            }
            catch (RuntimeException ex) {
                    System.out.println("Failed to print a pair");
                    throw ex;
            }
        }

        if (! on) System.out.println(" "+this +" is off");
        Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();

       while((pair = (PairOfHydrogenBondsElements)pairsIter.next()) != null){
             if ((PairOfHydrogenBondsElements.helixPair(pair.HOelement1) & PairOfHydrogenBondsElements.helixPair(pair.HOelement2)) |
                 (PairOfHydrogenBondsElements.sheetPair(pair.HOelement1) & PairOfHydrogenBondsElements.sheetPair(pair.HOelement2))){

                energyElement.set(pair);
                if(((HydrogenBondsPairsEnergyElement )energyElement).distanceAttributes1 != ((HydrogenBondsPairsEnergyElement )energyElement).HOelement1.getAttribute(HB_DistanceAttribute.key))
                    System.out.println("this is wird");
                if(energyElement .atoms().whereIs(atom) >= 0){
                    System.out.println(comment+" element-pair: ");
                    energyElement.atoms().print();

                //System.out.println(energyElement.oAtom1());
                    if(pair.lookAtThisPair){
                        //System.out.println(comment+": lookAtThisPair");
                        //System.out.println("getAtoms: "+ pair.getAtoms());
                        //System.out.println("size "+pair.getAtoms().size());
                        //System.out.println("pairs "+pair.getAtoms());
                    }
                    else{
                        System.out.println("dont look at this pair" +comment);
                        //pair.getAtoms().print();
                    }
                    //System.out.println("start Test pair !!!!!!!!!!!!!!!!!!");
                    try{
                        energyElement.test(totalEnergy,atom);
                    }
                    catch(SynchronizedUpdateException e)   {
                        System.out.println(e);

                    }
                    //System.out.println("End Test pair !!!!!!!!!!!!!!!!!!!!!");
                }
                ((HydrogenBondsPairsEnergyElement )energyElement).freeElenet();
             }
        }
        System.out.println("End Test "+comment);
    }
    /*
      OLD TEST 
    public void test(TotalEnergy totalEnergy, boolean testAll, boolean verbose) {
        double e;
        pairsOfHBEElementsList.print();
    	if (! reset) System.out.println(""+this +" is off");
    	double highestEnergy = -100000;
    	double lowestEnergy = 100000;
    	Iterator pairsIter = pairsOfHBEElementsList.iterator();
    	PairOfHydrogenBondsElements pair;
    	while ((pair  =  (PairOfHydrogenBondsElements)pairsIter.next()) != null) {
            if (!(DistanceMatrix.rMax()-pair.HOelement1.distance() >=0) || !(DistanceMatrix.rMax()-pair.HOelement2.distance() >=0)) continue;
            System.out.println("PairOfHydrogenBondsElements:"+pair);
            energyElement.setResidue(pair);
            e = energyElement.evaluate();
            if (e > highestEnergy) {
                System.out.println("********* Highest energy so far "+e+" ************");
                highestEnergy = e;
                energyElement.test(totalEnergy, verbose);	
            }
            else if(testAll) energyElement.test(totalEnergy, verbose);
            if (e < lowestEnergy) {
                System.out.println("********* Lowest energy so far "+e+" ************");
                lowestEnergy = e;
            }    		
    	}
    }
    */
}
