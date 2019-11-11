package meshi.energy.hydrogenBondsPairs;

import meshi.energy.hydrogenBond.HB_DistanceAttribute;
import meshi.energy.hydrogenBond.HydrogenBondsEnergyElement;
import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Utils;

public class HydrogenBondsPairsEnergyElement extends NonBondedEnergyElement {
    public static final double HELIX_FACTOR = 0.1;  // chen Sept. 2011
    public static final double PARALLEL_BETA_FACTOR = 1.8; //Chen 29.9.11 Trying to strengthen parallel beta that seems too week.
    static boolean toPrint = true;
    private final double width;
    private DistanceMatrix distanceMatrix;
    private PairOfHydrogenBondsElements.SsElement ssElement;
    //----------------------------------------- data ------------------------------------

    BetaParametersList betaParametersList;
    HelixParametersList helixParametersList;
    BetaParametersList betaBackgroundParametersList;
    HelixParametersList helixBackgroundParametersList;
    double punish,hpunish;
    private Atom oAtom1, hAtom1, oAtom2, hAtom2;
    private double weight;
    Distance HOelement1,HOelement2;
    HB_DistanceAttribute distanceAttributes1,distanceAttributes2;
    private double deDxOAtom1, deDyOAtom1, deDzOAtom1,deDxOAtom2, deDyOAtom2, deDzOAtom2;
    private double deDxHAtom1, deDyHAtom1, deDzHAtom1,deDxHAtom2, deDyHAtom2, deDzHAtom2;
    private double energy = 0;
    int ooDistance=-9999;
    int hhDistance = -9999;

    private boolean free = true;

    /*
    *  the probobility of a pair is::
    *  -log(P(D | T)) = - log { P{ T | D) * P (D) / P(T))
    *   where P(D) = the prior  prob to have 2 HB when picking 2 pairs of H-O getAtoms.Itdependes reset the number of residues that are in secondary strucutre segments
    *  yet, in thus version it will be a constant number
    *  given a specific pattern <x1,x2,x3,x4,x5,x6> P(T) is the number that this pattern was counted when we looked at two pairs of H-O getAtoms, in whicn at least
    *  one pair is bounded divided by the number of all patterns.
    *  given a specific pattern P(T | D) is the same as P(T), only now, both H-O pairs are bounded
    *
    * */


    private double PD_Beta=0.005, PD_Helix=0.0007; //P(D) for beta or helix
        private double PT_Beta_Background = 3987207, PT_Helix_Background = 57265603; //  the number of all patterns when we looked at two pairs of H-O getAtoms, in whicn at least  one pair is bounded
    private double PT_Beta_spesificPattern, PT_Helix_spesificPattern; //  given a specific pattern <x1,x2,x3,x4,x5,x6> P(T) is the number that this pattern was counted when we looked at two pairs of H-O getAtoms, in whicn at least one pair is bounded
    private double PT_givenD_Beta_Background = 54891,PT_givenD_Helix_Background=128659; //  the number of all patterns when we looked at two pairs of H-O getAtoms, in whicn both pairs are bounded
    private double PT_givenD_Beta_spesificPattern,PT_givenD_Helix_spesificPattern; //  given a specific pattern <x1,x2,x3,x4,x5,x6> P(T) is the number that this pattern was counted when we looked at two pairs of H-O getAtoms, in whicn both  pairs are bounded
    private double PT_Beta,PT_Helix; // PT_spesificPattern /  PT_Background
    private double PT_givenD_Beta, PT_givenD_Helix; //       PT_givenD_spesificPattern / PT_givenD_Background
    private double PD_givenT_Beta, PD_givenT_Helix;//    P(D | T) =  P(T | D) * P (D) / P(T)
    private double log_PD_givenT_Beta, log_PD_givenT_Helix;

    public double factor; //when factor >0 this patern has been seen in the datebase
                                                     //when factor <0 this patern is rare in the datebase
    private double pseodoCount =1;
    private int[] specialDisArray = null; //used when the user news one or more pairs of residues that share HB
    private boolean antiParallel = true;// in use only if specialDisArray is not null

    //----------------------------------------- constructors --------------------------------

    public HydrogenBondsPairsEnergyElement() {
        width = -1;
    }

    public HydrogenBondsPairsEnergyElement(HelixParametersList helixParametersList,
                                           BetaParametersList betaParametersList,
                                           HelixParametersList helixBackgroundParametersList,
                                           BetaParametersList betaBackgroundParametersList,
                                           double weight,
                                           double punish,
                                           double hpunish,
                                           double width,
                                           DistanceMatrix distanceMatrix) {
        this.helixParametersList= helixParametersList;
        this.betaParametersList = betaParametersList;
        this.helixBackgroundParametersList= helixBackgroundParametersList;
        this.betaBackgroundParametersList = betaBackgroundParametersList;
        this.weight = weight;
        this.punish = punish;
        this.hpunish = hpunish;
        this.width = width;
        this.distanceMatrix = distanceMatrix;
        Utils.println("HBP punihs: " + punish + " HBP hPunish: " + hpunish + "OOPunish = 100*punish   width=" + width);
    }

    public HydrogenBondsPairsEnergyElement(
                                           HelixParametersList helixParametersList,
                                           BetaParametersList betaParametersList,
                                           HelixParametersList helixBackgroundParametersList,
                                           BetaParametersList betaBackgroundParametersList,
                                           double weight,
                                           double punish,
                                           double hpunish,
                                           int[] specialDisArray,
                                           boolean antiParallel,
                                           double width,
                                           DistanceMatrix distanceMatrix) {
        this(helixParametersList, betaParametersList, helixBackgroundParametersList, betaBackgroundParametersList,weight, punish, hpunish, width,distanceMatrix);
        this.specialDisArray = specialDisArray ;
        this.antiParallel = antiParallel;
        if(specialDisArray != null)
            System.out.println("********************* Wornning: check that setResidue treat this case!!!");
    }


    //------------------------------------------ methods -------------------------------------

    public void set(Object obj){//obj should be PairOfHydrogenBondsElements
        free = false;
        PairOfHydrogenBondsElements pair = (PairOfHydrogenBondsElements) obj;
        if(! pair.lookAtThisPair)
            throw new RuntimeException("problem in HydrogenBondsPairsEnergyElement:this pair should not be looked at: "+pair);

        atoms = pair.atoms();
        this.HOelement1 = pair.HOelement1;
        this.HOelement2 = pair.HOelement2;
        distanceAttributes1 = (HB_DistanceAttribute)(HOelement1.getAttribute(HB_DistanceAttribute.key));
        distanceAttributes2 = (HB_DistanceAttribute)(HOelement2.getAttribute(HB_DistanceAttribute.key));
        if (distanceAttributes1 == null || distanceAttributes2 == null)
                   throw new RuntimeException("problem in HydrogenBondsPairsEnergyElement: one of the element is null");
        oAtom1 = distanceAttributes1.oAtom;
        hAtom1 = distanceAttributes1.hAtom;
        oAtom2 = distanceAttributes2.oAtom;
        hAtom2 = distanceAttributes2.hAtom;

        hhDistance = pair.hhSeqDistance() ;//parameters.h1h2SeqGap;
        ooDistance = pair.ooSeqDistance() ;//parameters .o1o2SeqGap;

       // if (hhDistance == 0){
        //    factor = hpunish; //TODO seperate between beta and helix ??
       // }
        //TODO special tretment for the case of ooDistance == 0 ??
      //  else {
                double parameterValue, backgroundParameterValue;
                HydrogenBondsPairsParameters parameters;
                HydrogenBondsPairsParameters backgroundParameters;


                if(pair.getPairValue() == -1000 | pair.getPairBackgroundValue() == -1000){
                            if(pair.isHelixPair() ) {
                                     //System.out.println("HydrogenBondsPairsEnergyElement: helixPair "+pair);

                                    ssElement = PairOfHydrogenBondsElements.SsElement.HELIX;
                                    parameters = (HydrogenBondsPairsParameters)helixParametersList.parameters(pair);
                                     backgroundParameters =  (HydrogenBondsPairsParameters)helixBackgroundParametersList.parameters(pair);
                                    parameterValue = parameters.value();
                                    backgroundParameterValue = backgroundParameters.value();
                                    pair.setPairValue(parameterValue );
                                    pair.setPairBackgroundValue(backgroundParameterValue);
                                    if(parameterValue == HydrogenBondsPairsParametersList.nonExsistValue || backgroundParameterValue == HydrogenBondsPairsParametersList.nonExsistValue ){
                                          factor = punish;
                                          ssElement = PairOfHydrogenBondsElements.SsElement.HELIX_NONEXISTENT;
                                      }
                                    else{
                                    PT_Helix_spesificPattern = backgroundParameterValue;
                                    PT_Helix = PT_Helix_spesificPattern / PT_Helix_Background;
                                    PT_givenD_Helix_spesificPattern = pseodoCount+parameterValue;
                                    PT_givenD_Helix = PT_givenD_Helix_spesificPattern / PT_givenD_Helix_Background;
                                    PD_givenT_Helix = PT_givenD_Helix  / PT_Helix;  // PT_givenD_Helix * PD_Helix / PT_Helix;
                                    log_PD_givenT_Helix = Math.log(PD_givenT_Helix);
                                    factor = -1*log_PD_givenT_Helix*HELIX_FACTOR; //HELIX_FACTOR was added by Chen, Sept. 2011
                                    }
                                }
                                 else if(pair.isBetaPair() ) {
                                      ssElement = PairOfHydrogenBondsElements.SsElement.BETA;
                                      parameters = (HydrogenBondsPairsParameters)betaParametersList.parameters(pair ) ;
                                      backgroundParameters =  (HydrogenBondsPairsParameters)betaBackgroundParametersList.parameters(pair);
                                      parameterValue = parameters.value();
                                      backgroundParameterValue = backgroundParameters.value();
                                      pair.setPairValue(parameterValue );
                                      pair.setPairBackgroundValue(backgroundParameterValue);
                                      if(parameterValue == HydrogenBondsPairsParametersList.nonExsistValue || backgroundParameterValue == HydrogenBondsPairsParametersList.nonExsistValue ){
                                          factor = punish;
                                          ssElement = PairOfHydrogenBondsElements.SsElement.BETA_NONEXISTENT;
                                      }
                                      else {
                                              PT_Beta_spesificPattern = backgroundParameterValue;
                                              PT_Beta = PT_Beta_spesificPattern / PT_Beta_Background;
                                              PT_givenD_Beta_spesificPattern = pseodoCount+parameterValue;
                                              PT_givenD_Beta = PT_givenD_Beta_spesificPattern / PT_givenD_Beta_Background;
                                              PD_givenT_Beta = PT_givenD_Beta / PT_Beta;  // PT_givenD_Beta * PD_Beta / PT_Beta;
                                              log_PD_givenT_Beta = Math.log(PD_givenT_Beta);
                                              factor = -1*log_PD_givenT_Beta;
                                              if (factor < 0) {
                                                  if (parameters.isParallel())  {
                                                    ssElement = PairOfHydrogenBondsElements.SsElement.PARALLEL_BETA;
                                                    factor *= PARALLEL_BETA_FACTOR;
                                                  }
                                                  else ssElement = PairOfHydrogenBondsElements.SsElement.ANTI_PARALLEL_BETA;
                                              }
                                          }
                                      }

                                 else if(pair. isOneBetta() )
                                      throw new RuntimeException("HydrogenBondsEnergyElement: oneBeta");
                                 else if(pair.isOneHelix() )
                                      throw new RuntimeException("HydrogenBondsEnergyElement: oneHelix");
                                 else
                                       throw new RuntimeException("pair should be beta/helix pair\n"+pair);
                            pair.setPairFactor(factor, ssElement);
                }
                else
                    factor = pair.getPairFactor();
            }


      //  }






        /*
        if (pair.seqDistance1() == -4 & //parameters.firstHBSeqGap == -4 &
                pair.seqDistance2() == -4 & //parameters .secondHBSeqGap == -4 &
                hhDistance == 1&//parameters .h1h2SeqGap == 1 &
                ooDistance == 1 &//parameters .o1o2SeqGap == 1 &
                pair.hoSeqDistance() == 5&//parameters .h1o2SeqGap == 5 &
                pair.ohSeqDistance() ==-3)//parameters .o1h2SeqGap == -3)  /
                factor = 10; //This is a helix !
        else if  (pair.isHelixPair() &
                pair.seqDistance1() == -4 & //parameters.firstHBSeqGap == -4 &
                pair.seqDistance2() == -4 & //parameters .secondHBSeqGap == -4 &
                hhDistance ==  ooDistance )
                factor = 1;  //This is a helix !       changed reset 27/11/06 21:00

                 else if(ooDistance ==  0) {
                     //if(Math.abs((pair. hAtom1.residueNumber() - pair. hAtom2.residueNumber())) > 5)                  //TODO add a spacial punish to this case, if needed (if results indicates that there are more then expected o getAtoms with more the one hbonds
                           // System.out.println("HBPEnergyElement: same o atom: "+pair);
                     //System.out.println("punish "+punish);
                  factor = 100 * punish+1;
                  }

            else if (hhDistance == 0){
                if(pair.sheetPair(HOelement1 ) | pair .sheetPair(HOelement2 ) )
                      factor = 100*punish+1;//hpunish;//*10;
                else factor = punish;//hpunish;          //negative
            }
            else  {

                double parameterValue ;
                HydrogenBondsPairsParameters parameters = null;
                if(pair.getPairValue() != -1000)
                                      parameterValue = pair.getPairValue();
                else{
                               if(pair.isHelixPair() ) {
                                     System.out.println("HydrogenBondsPairsEnergyElement: helixPair "+pair);
                                     parameters = (HydrogenBondsPairsParameters)helixParametersList.parameters(pair);
                                }
                                 else if(pair.isBetaPair() )
                                      parameters = (HydrogenBondsPairsParameters)betaParametersList.parameters(pair ) ;
                                 else if(pair. isOneBetta() )
                                      System.out.println("HydrogenBondsEnergyElement: oneBeta");
                                 else if(pair.isOneHelix() )
                                      System.out.println("HydrogenBondsEnergyElement: oneHelix");
                                 else
                                       throw new RuntimeException("pair should be beta/helix pair");
                                 parameterValue = parameters.value();
                                 pair.setPairValue(parameterValue );
                 }


         if (parameterValue < 0){
            if (ooDistance ==  0)
                      System.out .println("ooDistance ==  0 :"+pair);
            if(pair.isHelixPair())
                factor = punish*10;//negative
            else
                factor = punish; //negative
        }
         else if(specialDisArray != null){
            if(antiParallel){
                if (!antiParallel(parameters ))
                    factor = 0;
                else if(antiParallelshift(hAtom1 .residueNumber(),oAtom1 .residueNumber()))
                    factor = -1;
                else
                    factor = parameterValue/10;
            }
            else
                if( ! parallel(parameters ) )
                    factor = 0;
                else if(parallelshift(hAtom1 .residueNumber(),oAtom1 .residueNumber()))
                    factor = -1;
                else
                    factor = parameterValue/10;
        }
         else
            factor = parameterValue/10;//positive
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");

      }  */
       //TODO treat cerfully the specialDisArry
      /*   This isthe old trearment!!
        if(specialDisArray != null){
            if(antiParallel){
                if (!antiParallel(parameters ))
                    factor = 0;
                else if(antiParallelshift(hAtom1 .residueNumber(),oAtom1 .residueNumber()))
                    factor = -1;
                else
                    factor = parameterValue/10;
            }
            else
                if( ! parallel(parameters ) )
                    factor = 0;
                else if(parallelshift(hAtom1 .residueNumber(),oAtom1 .residueNumber()))
                    factor = -1;
                else
                    factor = parameterValue/10;
          */




    private boolean antiParallel(HydrogenBondsPairsParameters parameters){
        return (parameters .h1h2SeqGap ==-1*parameters.o1o2SeqGap & parameters.h1o2SeqGap == -1*parameters.o1h2SeqGap) ;
    }

    private boolean parallel(HydrogenBondsPairsParameters parameters){
        return  ((parameters .h1h2SeqGap == parameters.o1o2SeqGap & (parameters .h1h2SeqGap == 2 | parameters .h1h2SeqGap == 4) & parameters.h1o2SeqGap == -1*parameters.o1h2SeqGap) |
                          (parameters .h1h2SeqGap == -1*parameters.o1o2SeqGap & parameters.h1o2SeqGap == parameters.o1h2SeqGap +2  & (parameters.h1o2SeqGap == 4 | parameters.h1o2SeqGap == 2 | parameters.h1o2SeqGap == 0 | parameters.h1o2SeqGap == -2)));

    }

    private boolean antiParallelshift (int resNum1,int resNum2){
        return Math.abs(resNum1 - resNum2) % 2 != Math.abs(specialDisArray[0] - specialDisArray[1]) % 2;
    }

     private boolean parallelshift (int resNum1,int resNum2){
        return Math.abs(resNum1 - resNum2) % 2 != Math.abs(specialDisArray[0] - (specialDisArray[1]+1)) % 2;
    }

    public void freeElenet(){
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        free = true;
        atoms = null;
        HOelement1 = HOelement2 = null;
        distanceAttributes2 = distanceAttributes1 = null;
        factor = 1/0.0;
        deDxHAtom1 = deDxHAtom2 =deDxOAtom1 =deDxOAtom2 = 1/0.0;
        deDyHAtom1 = deDyHAtom2 =deDyOAtom1 =deDyOAtom2 =1/0.0;
        deDzHAtom1 =deDzHAtom2 =deDzOAtom1 =deDzOAtom2 =1/0.0;
        energy =1/0.0;
        oAtom1 =oAtom2 =hAtom1 =hAtom2=null;
    }

    /* (non-Javadoc)
    * @see meshi.energy.EnergyElement#evaluate()
    */
    public double evaluate() {
       // int   x;
      //  if(hAtom1.number() ==29 | hAtom2.number() == 29)
      //      x = 1;
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        updateEnergy();
        updateAtoms();
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        return energy * weight;
    }

    public double updateEnergy() {
        HydrogenBondsEnergyElement hbElement = new HydrogenBondsEnergyElement(1, width,distanceMatrix);     //TODO change this
        hbElement.set(HOelement1);
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        if(distanceAttributes1 != hbElement.hb_Attribute )
            System.out.println("this is wird");
        hbElement.updateEnergy();
         hbElement.freeElement();
        hbElement.set(HOelement2); //TODO try use new element
         hbElement.updateEnergy();
         hbElement.freeElement();
        //TODO take care of this in more elegant way, now there are double calculations of HB energy

        double e1, e2;
        double dE1DxOAtom1, dE1DyOAtom1, dE1DzOAtom1;//, dE1DxOAtom2, dE1DyOAtom2, dE1DzOAtom2;
        double dE2DxOAtom2, dE2DyOAtom2, dE2DzOAtom2;//, dE2DxOAtom1, dE2DyOAtom1, dE2DzOAtom1;
        double dE1DxHAtom1, dE1DyHAtom1, dE1DzHAtom1;//, dE1DxHAtom2, dE1DyHAtom2, dE1DzHAtom2;
        double dE2DxHAtom2, dE2DyHAtom2, dE2DzHAtom2;//, dE2DxHAtom1, dE2DyHAtom1, dE2DzHAtom1;

        e1 = distanceAttributes1.energy;

        dE1DxOAtom1 = distanceAttributes1.deDxOAtom;
        dE1DyOAtom1 = distanceAttributes1.deDyOAtom;
        dE1DzOAtom1 = distanceAttributes1.deDzOAtom;
        dE1DxHAtom1 = distanceAttributes1.deDxHAtom;
        dE1DyHAtom1 = distanceAttributes1.deDyHAtom;
        dE1DzHAtom1 = distanceAttributes1.deDzHAtom;

        e2 = distanceAttributes2.energy;

        dE2DxOAtom2 = distanceAttributes2.deDxOAtom;
        dE2DyOAtom2 = distanceAttributes2.deDyOAtom;
        dE2DzOAtom2 = distanceAttributes2.deDzOAtom;
        dE2DxHAtom2 = distanceAttributes2.deDxHAtom;
        dE2DyHAtom2 = distanceAttributes2.deDyHAtom;
        dE2DzHAtom2 = distanceAttributes2.deDzHAtom;
        deDxOAtom1 =  factor * (e2 * dE1DxOAtom1);  //*-1 22 Oct 2010 (Ami)
        deDyOAtom1 =  factor * (e2 * dE1DyOAtom1);
        deDzOAtom1 =  factor * (e2 * dE1DzOAtom1);
        deDxHAtom1 =  factor * (e2 * dE1DxHAtom1);
        deDyHAtom1 =  factor * (e2 * dE1DyHAtom1);
        deDzHAtom1 =  factor * (e2 * dE1DzHAtom1);

        deDxOAtom2 =  factor * (e1 * dE2DxOAtom2 );     //*-1 22 Oct 2010 (Ami)
        deDyOAtom2 =  factor * (e1 * dE2DyOAtom2 );
        deDzOAtom2 =  factor * (e1 * dE2DzOAtom2 );
        deDxHAtom2 =  factor * (e1 * dE2DxHAtom2 );
        deDyHAtom2 =  factor * (e1 * dE2DyHAtom2 );
        deDzHAtom2 =  factor * (e1 * dE2DzHAtom2 );

        // e1, e2 are negative values => negative (changed from positive, reset 22 Oct, 2010, by Ami) factor means negative (good) energy
        //                               positive (changed from negative, reset 22 Oct, 2010, by Ami) factor means positive (bad) energy
        energy = factor * e1 * e2;              //*-1 22 Oct 2010 (Ami)
        if (factor > 0 && energy <0)  // factor < 0 changed reset 22 Oct 2010 (Ami)
            System.out.println("HBPEnergy Element: e1*e2 is negative !! energy: "+energy+" e1:  "+e1+" e2: "+e2+" weight: "+weight+" factor: "+factor);
        if (energy < -10000 || energy > 10000){
            System.out.println("HBPEnergy Element: energy is very low or vey high: energy: "+energy+" e1:  "+e1+" e2: "+e2+" weight: "+weight+" factor: "+factor);
            throw new RuntimeException(""+this);
        }
        return energy;
    }

    public void updateAtoms(){
        if (! oAtom1.frozen()) {
            oAtom1.addToFx(-1 * deDxOAtom1 * weight); // force = -derivative
            oAtom1.addToFy(-1 * deDyOAtom1 * weight); // force = -derivative
            oAtom1.addToFz(-1 * deDzOAtom1 * weight); // force = -derivative       
        }
        if (! hAtom1.frozen()) {
            hAtom1.addToFx(-1 * deDxHAtom1 * weight);
            hAtom1.addToFy(-1 * deDyHAtom1 * weight);
            hAtom1.addToFz(-1 * deDzHAtom1 * weight);
        }
        if (! oAtom2.frozen()) {
            oAtom2.addToFx(-1 * deDxOAtom2 * weight); // force = -derivative
            oAtom2.addToFy(-1 * deDyOAtom2 * weight); // force = -derivative
            oAtom2.addToFz(-1 * deDzOAtom2 * weight); // force = -derivative            
        }
        if (! hAtom2.frozen()) {
            hAtom2.addToFx(-1 * deDxHAtom2 * weight);
            hAtom2.addToFy(-1 * deDyHAtom2 * weight);
            hAtom2.addToFz(-1 * deDzHAtom2 * weight);
        }
    }

    /* (non-Javadoc)
    * @see meshi.energy.EnergyElement#setAtoms()
    */
    protected void setAtoms() {
        throw new UnsupportedOperationException("setAtoms() may not be used by HydrogenBondsEnergyElement for "+
                                   "efficiency.") ;
    }

    public String toString() {
        if(free) return "HydrogenBondsPairsEnergyElement: is free (no values)";
        return "HydrogenBondsPairsEnergyElement: oAtom1 "+oAtom1.residueNumber()+" HAtom1 " +hAtom1.residueNumber()+" oAtom2 "+oAtom2.residueNumber()+" HAtom2 " +hAtom2.residueNumber();
        //return "HydrogenBondsPairsEnergyElement: factor "+factor+" "+pair;
    }
}
