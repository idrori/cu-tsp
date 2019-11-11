package meshi.energy.hydrogenBond;

import meshi.util.MeshiAttribute;
import meshi.molecularElements.atoms.*;

public class HB_AtomAttribute implements MeshiAttribute{
    public final boolean isO;
    public final boolean isH;
    public static final int key = HYDROGEN_BONDS_ATTRIBUTE;
    public HB_AtomAttribute(boolean isH,boolean isO) {
	if ((isO & isH) | (!((isO | isH)))) 
	    throw new RuntimeException("Weird parameters "+isO+" "+isH); 
	this.isH = isH;
	this.isO = isO;
    }
    public final int key(){return key;}
    public String toString() {
	    	return "HB_AtomAttribute "+key+" "+isO+" "+isH;
    }
}
