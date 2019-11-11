/*
 * Created on 16/11/2004
 * Window - Preferences - Java - Code Style - Code Templates
 */
package meshi.energy.hydrogenBond;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.atoms.*;
import meshi.geometry.Distance;

/**
 * @author amilev
 *
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class HydrogenBondsParametersList extends LennardJonesParametersList {
	AtomList nitrogens;
	
	public HydrogenBondsParametersList(String parametersFileName,AtomList nitrogens){
		super(parametersFileName);
		this.nitrogens = nitrogens;
	}

    public HydrogenBondsParametersList(String parametersFileName){
		super(parametersFileName);
    }
	
	/* (non-Javadoc)
	 * @see meshi.energy.ParametersList#getParameters(java.lang.Object)
	 */
	public Parameters parameters(Object obj) {
		Distance distance = (Distance) obj;
		int largeType = distance.largeType.ordinal();
		int smallType = distance.smallType.ordinal();
		try {
			return (HydrogenBondsParameters) get(largeType*(largeType-1)/2+smallType);
	}	
		catch (Exception e) {
			throw new RuntimeException("Cannot find parameters for:\n"+
						   distance+"\n"+
						   distance.atom1()+"\n"+
						   distance.atom2()+"\n"+
						   largeType+"-"+distance.largeType+" "+ 
						   smallType+"-"+distance.smallType+" "+size()+" "+
						   (largeType*(largeType-1)/2+smallType)+"\n"+e); 
		}	
	}

	/* (non-Javadoc)
	 * @see meshi.energy.ParametersList#createParameters(java.lang.String)
	 */
	public Parameters createParameters(String line) {
			return new HydrogenBondsParameters(new StringTokenizer(line));
	}

	/**
	 * @return nitrogens
	 */
	public AtomList nitrogens() {
		return nitrogens;
	}

}
