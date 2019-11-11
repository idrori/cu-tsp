package meshi.energy.hydrogenBond;

import meshi.util.MeshiAttribute;
import meshi.molecularElements.atoms.*;

//--------------------------------------------------------------------------------------------
    
    //------------------------------- class DistanceAttribute -----------------------------

    //--------------------------------------------------------------------------------------------

public class HB_DistanceAttribute implements MeshiAttribute{
    public double energy;	
    public double deDxOAtom;
    public double deDyOAtom;
    public double deDzOAtom;
    public double deDxHAtom;
    public double deDyHAtom;
    public double deDzHAtom;
    public Atom oAtom, hAtom;
    public final boolean isHbond;
    public double epsilon, sigma, sigma6, sigma6EpsilonFour, minusTwelveSigma6;
    public HydrogenBondsParameters parameters;
    public static final int key = HYDROGEN_BONDS_ATTRIBUTE;

    public HB_DistanceAttribute (boolean isHbond){
	this.isHbond = isHbond;
    }

    public void setParameters(HydrogenBondsParameters parameters){
        epsilon = parameters.epsilon;
        sigma = parameters.sigma;
        sigma6 = parameters.sigma6;
        sigma6EpsilonFour = parameters.sigma6EpsilonFour;
        minusTwelveSigma6 = parameters.minusTwelveSigma6;
    }
    
    public void set(double deDxOAtom,double deDyOAtom,double deDzOAtom,
			double deDxHAtom,double deDyHAtom,double deDzHAtom,
			double energy,Atom oAtom,Atom hAtom){
	set(deDxOAtom,  deDyOAtom, deDzOAtom,
			deDxHAtom, deDyHAtom, deDzHAtom,
			 energy);
    set(oAtom, hAtom);   
    }

    public void set(double deDxOAtom,double deDyOAtom,double deDzOAtom,
			double deDxHAtom,double deDyHAtom,double deDzHAtom,
			double energy){
	this.deDxOAtom = deDxOAtom;
	this.deDyOAtom = deDyOAtom;
	this.deDzOAtom = deDzOAtom;
	this.deDxHAtom = deDxHAtom;
	this.deDyHAtom = deDyHAtom;
	this.deDzHAtom = deDzHAtom;
	this.energy = energy;
    }

    public void set(Atom atom1,Atom atom2){//one should be oAtom and the other hAtom
	HB_AtomAttribute atom11Attribute =
	    (HB_AtomAttribute) atom1.getAttribute(HB_AtomAttribute.key);
	if (atom11Attribute == null) {
	    System.out.println("atom11Attribute is null \n"+atom1+" "+ atom1.getAttribute(HB_AtomAttribute.key)+"\n"+
			       atom2+" "+ atom2.getAttribute(HB_AtomAttribute.key)); 
	    throw new RuntimeException("XXXX");
	}
        if (atom11Attribute.isH) {
            hAtom = atom1;
            oAtom = atom2;
        }
        else {
            hAtom = atom2;
            oAtom = atom1;
        }
    }


    
    
    public final int key() {return key;}
}
