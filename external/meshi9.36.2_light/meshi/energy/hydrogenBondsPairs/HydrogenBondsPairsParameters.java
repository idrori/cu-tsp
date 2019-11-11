/*
 * Created on 26/12/2004
 */
package meshi.energy.hydrogenBondsPairs;

import meshi.energy.Parameters;
import meshi.util.*;
import java.util.StringTokenizer;

/**
 * @author amilev
 *
 * Let di be the residue number of a hydrogen-donor at hydrogen-bond i,
 *  and ai the residue number of a hydrogen-acceptors in the same hydrogen-bond.
 *  Thus, we can uniquely represent each pattern <hbi ,hbk> by sextuplet:  <|ai-di| ,|ak-dk|,|di-dk|,  |di-ak|,|ai-dk|,|ai-ak|>.
 *  To reduce the number of possible sextuplets, |ai-di| is bounded by 6, and|ai-dk| is bounded by 10.

 **/
public class HydrogenBondsPairsParameters implements Comparable, Parameters  {

    //Chen 29.9.11 trying to strengthen parallel beta, that apparently our statistics does not like
    public static final HydrogenBondsPairsParameters parallelBeta1 = new HydrogenBondsPairsParameters(11,-11,6,-6,2,0);
    public static final HydrogenBondsPairsParameters parallelBeta2 = new HydrogenBondsPairsParameters(11,-11,6,-6,0,-2);
    public static final HydrogenBondsPairsParameters parallelBeta3 = new HydrogenBondsPairsParameters(-11,-11,2,2,6,-6);
    public static final HydrogenBondsPairsParameters parallelBeta4 = new HydrogenBondsPairsParameters(11,11,2,2,-6,6);
    public static final HydrogenBondsPairsParameters parallelBeta5 = new HydrogenBondsPairsParameters(11,-11,6,-6,4,2);
    public static final HydrogenBondsPairsParameters parallelBeta6 = new HydrogenBondsPairsParameters(11,-11,6,-6,-2,-4);
    public static final HydrogenBondsPairsParameters parallelBeta7 = new HydrogenBondsPairsParameters(11,11,4,4,-6,6);

    //------------------------------------ data ------------------------------------
    /**
     * Key values:
     * firstHBSeqGap = the gap between the two getAtoms in the first HB   (O1-H1)
     * secondHBSeqGap = the gap between the two getAtoms in the second HB   (O2-H2)
     * h1h2SeqGap = the gap between the two Hydrogens
     * o1o2SeqGap = the gap between the two Oxygens
     * h1o2SeqGap = the gap between the H of the first HB and the O of the second HB
     * o1h2SeqGap = the gap between the O of the first HB and the H of the second HB
     * value is proportionaly to the friquency of such pair of HBs in the data base
     **/
	public  int firstHBSeqGap=0,secondHBSeqGap=0,h1h2SeqGap=0,o1o2SeqGap=0,h1o2SeqGap=0,o1h2SeqGap=0;

    private int value=0;
    public final int value(){return value;}
    //------------------------------------ constructors ----------------------------

    /*
     * create an element from line that contains all the key values and the friquncy value
     */
    public HydrogenBondsPairsParameters(String line) {

            this(line, new StringTokenizer(line,": "));



        }
	    
    private HydrogenBondsPairsParameters(String wholeLine, StringTokenizer line) {
        if(line.countTokens() == 14 )            {
       //     throw new RuntimeException("There is a problem in the helix  parameters file; \n " +
        //            "each row should look like: \" firstBond : X   secondBond : X    h_dis : X    o_dis : X    ho-dis : X    oh_dis : X *** X \" \n" +
         //           "but this line has a problem: "+wholeLine);
        line.nextToken();
        firstHBSeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        secondHBSeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        h1h2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        o1o2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        h1o2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        o1h2SeqGap = Utils.toInt(line.nextToken());
        line.nextToken();
        value = Utils.toInt(line.nextToken());
            }
    }
	    
    public HydrogenBondsPairsParameters(int firstHBGap,int secondHBGap,int h1h2Gap,int o1o2Gap,int h1o2Gap,int o1h2Gap){
        firstHBSeqGap =firstHBGap;
        secondHBSeqGap = secondHBGap;
        h1h2SeqGap = h1h2Gap;
        o1o2SeqGap = o1o2Gap;
        h1o2SeqGap = h1o2Gap;
        o1h2SeqGap = o1h2Gap;
    }

    public HydrogenBondsPairsParameters(int firstHBGap,int secondHBGap,int h1h2Gap,int o1o2Gap,int h1o2Gap,int o1h2Gap,int value){
        this(firstHBGap,secondHBGap,h1h2Gap,o1o2Gap,h1o2Gap,o1h2Gap);
        this.value = value;
    }

    //--------------------------------- methods ------------------------------------------


    public boolean isParallel() {
       return equals(parallelBeta1)|| equals(parallelBeta2)|| equals(parallelBeta3) || equals(parallelBeta4)||
              equals(parallelBeta5)|| equals(parallelBeta5)|| equals(parallelBeta6)|| equals(parallelBeta7);
    }
    public boolean equals(Object obj){
      HydrogenBondsPairsParameters other = (HydrogenBondsPairsParameters) obj;
        if (firstHBSeqGap != other.firstHBSeqGap)return false;
        if (secondHBSeqGap != other.secondHBSeqGap) return false;
        if (h1h2SeqGap != other.h1h2SeqGap) return false;
        if (o1o2SeqGap != other.o1o2SeqGap) return false;
        if (h1o2SeqGap != other.h1o2SeqGap) return false;
        if (o1h2SeqGap != other.o1h2SeqGap) return false;
        return true;
    }
    public String toString()
    {
        return "firstBond: "+firstHBSeqGap+"   secondBond: "+secondHBSeqGap+"    h_dis: "+h1h2SeqGap+"    o_dis: "+o1o2SeqGap+"    ho-dis: "+h1o2SeqGap+"    oh_dis: "+o1h2SeqGap+" value: "+value;
    }
	    
	   
    public int compareTo(Object obj) {
        HydrogenBondsPairsParameters other = (HydrogenBondsPairsParameters)obj;
		if (firstHBSeqGap == other.firstHBSeqGap) {
		    if (secondHBSeqGap == other.secondHBSeqGap) {
                if (h1h2SeqGap == other.h1h2SeqGap) {
                    if (o1o2SeqGap == other.o1o2SeqGap) {
                        if (h1o2SeqGap == other.h1o2SeqGap) {
                            if (o1h2SeqGap == other.o1h2SeqGap) return 0;
                            else {
                                if (o1h2SeqGap > other.o1h2SeqGap) return 1;
                                return -1;
                            }
                        }
                        else {
                            if (h1o2SeqGap > other.h1o2SeqGap) return 1;
                            return -1;
                        }
                    }
                    else {
                        if (o1o2SeqGap > other.o1o2SeqGap) return 1;
                        return -1;
                    }
                }
                else {
                    if (h1h2SeqGap > other.h1h2SeqGap) return 1;
                    return -1;
                }
		    }
		    else {
                if (secondHBSeqGap > other.secondHBSeqGap) return 1;
                return -1;
		    }
		}
		else {
		    if (firstHBSeqGap > other.firstHBSeqGap) return 1;
		    return -1;
		}
    }


	    
}
