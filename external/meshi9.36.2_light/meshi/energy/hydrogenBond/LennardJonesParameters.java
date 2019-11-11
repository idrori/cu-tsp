package meshi.energy.hydrogenBond;
import meshi.energy.*;
import meshi.util.*;
import meshi.util.filters.*;
import java.util.*;

public class LennardJonesParameters implements MeshiAttribute, Parameters {
    public final double epsilon, sigma, sigma6, sigma6EpsilonFour, minusTwelveSigma6;  
    public static final int key = LENNARD_JONES_ELEMENT_ATTRIBUTE;
    public LennardJonesParameters() {
	epsilon = sigma = sigma6 = sigma6EpsilonFour = minusTwelveSigma6 = -1;
    }

    public LennardJonesParameters(StringTokenizer st) {
	epsilon = Utils.toDouble(st.nextToken());
	sigma = Utils.toDouble(st.nextToken());
	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
	sigma6EpsilonFour = 4*epsilon*sigma6;
	minusTwelveSigma6 = -12*sigma6;
    }
    
    public Filter isA() {return (new isA());}
    public final int key() {return key;}
    public Parameters create(StringTokenizer stringTokenizer) {
	return (new  LennardJonesParameters(stringTokenizer));
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof LennardJonesParameters);
	}
    }

    public String toString() {return ""+epsilon+"\t"+sigma;}
}
