package meshi.energy.hydrogenBondsPairs;

import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.ParametersList;

public abstract class HydrogenBondsPairsParametersList  extends ParametersList{

   public static final int nonExsistValue = -5;
    private final Parameters nonExsisting = new HydrogenBondsPairsParameters(1,1,1,1,1,1,nonExsistValue);

    public int searchCounter =0;

    //---------------------------------- constructor -------------------------------
    public HydrogenBondsPairsParametersList (){
        super();                //create an empty list with filter IsParameter
    }



    public HydrogenBondsPairsParametersList(String parametersFileName) {
        super(parametersFileName ,true);  //true means it is sortable parameter list using compare
    }

    //---------------------------------- methods -----------------------------------

     public Parameters parameters(Object baseElement) {
        searchCounter ++;
        Parameters key = getKey(baseElement);
        Parameters params =   getParameters(key);
         //Chen 29.9.11 changed as part of my effort to strengthen parallel beta
         if (params == null)
             return nonExsisting ;
         return params;
     }

    /**
     *
     * @param baseElement should be instance of PairOfHydrogenBondsElements
     * @return HydrogenBondsPairsParameters with 6 values to use un serching the freqency value
     */
    public Parameters getKey(Object baseElement) {
        PairOfHydrogenBondsElements element = (PairOfHydrogenBondsElements)baseElement;
        return new HydrogenBondsPairsParameters(element.seqDistance1,
                                                element.seqDistance2,
                                                element.hhSeqDistance,
                                                element.ooSeqDistance,
                                                element.hoSeqDistance,
                                                element.ohSeqDistance);
    }

    public Parameters createParameters(String line) {
        return new HydrogenBondsPairsParameters(line);
    }

 }


