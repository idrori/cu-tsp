
package meshi.energy.hydrogenBond;
import meshi.util.*;
import meshi.molecularElements.*;   
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import java.util.*;



public class LennardJonesParametersList  extends ParametersList{
        public LennardJonesParametersList(String parametersFileName) {
	    super(parametersFileName,false);
	}

     public Parameters createParameters(String line) {
	return new LennardJonesParameters(new StringTokenizer(line));
    }
    
    public Parameters parameters(Object obj) {
	Distance distance = (Distance) obj;
	int largeType = distance.largeType.ordinal();
	int smallType = distance.smallType.ordinal();
	try {
	    return (LennardJonesParameters) get(largeType*(largeType+1)/2+smallType);
	}
	catch (Exception e) {
	    throw new RuntimeException("Error in  LennardJonesParametersList.parameters()\n"+
				       "Failed to find parameters for "+obj+"\n"+
				       "large type = "+largeType+" smallType = "+ 
				       smallType+"\n"+e); 
	}	
   }
}
