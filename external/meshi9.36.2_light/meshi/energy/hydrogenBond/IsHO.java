/*
 * Created on 20/12/2004
 *
 * Window - Preferences - Java - Code Style - Code Templates
 */
package meshi.energy.hydrogenBond;

import java.util.Arrays;

import meshi.molecularElements.atoms.*;
import meshi.geometry.Distance;
import meshi.parameters.*;
import meshi.util.filters.Filter;

/**
 * @author amilev
 *
 * Window - Preferences - Java - Code Style - Code Templates
 */

public class IsHO implements Filter {
    public boolean accept(Object obj){
        Distance dis = (Distance)obj;
        HB_AtomAttribute atom1Attribute =
            (HB_AtomAttribute) (dis.atom1().getAttribute(HB_AtomAttribute.key));
    	HB_AtomAttribute atom2Attribute =
            (HB_AtomAttribute) (dis.atom2().getAttribute(HB_AtomAttribute.key));
	if ( atom1Attribute == null || atom2Attribute == null)//meens that at list one of the atom is not H or O
	    //since we add attribute just to H or O
	    return false;
        return ( ( atom1Attribute.isH &&
                   atom2Attribute.isO )
                 ||  
                 ( atom1Attribute.isO &&
                   atom2Attribute.isH ) );
        		
    }

}
