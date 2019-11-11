package meshi.energy.rapdf;
import meshi.energy.EnergyElement;
import meshi.energy.TotalEnergy;
import meshi.geometry.Coordinates;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;


public class RapdfEnergyElement extends EnergyElement {
    protected Atom atom, evaluatedAtom;   
    protected double force, forceSave;
    protected double evaluatedX;
    protected double evaluatedY;
    protected double evaluatedZ;
    //protected boolean frozen;
    
    protected FreeDistance distance;
    
    public RapdfEnergyElement() {}
    
    public RapdfEnergyElement(Atom inputAtom,double x,double y,double z, double weight) {
        evaluatedX=x;
        evaluatedY=y;
        evaluatedZ=z;
	
	atom = inputAtom;
	atoms = new AtomList(inputAtom.molecularSystem);
	atoms.add(atom);
        reset(evaluatedX,evaluatedY,evaluatedZ);
	force = forceSave = weight;
	frozen = atom.frozen();
    }


    /*** do nothing in the test stage,
     this function is not ment to compute derivatives for ***/
    public void test(TotalEnergy totalEnergy, Atom atom)
        {
            return;
        }

    protected void reset() {
        reset(atom.x(),atom.y(),atom.z());
    }
    protected void reset(double x, double y, double z) {
        evaluatedAtom =  new Atom("eval", null, AtomType.XXX, new Coordinates(x,y,z),1, 0.0,new MolecularSystem());

        distance = new FreeDistance(atom,evaluatedAtom);
    }
    protected void setAtoms(){};

    public void remove() {
        force = 0;
    }

    /*** compute rapdf score of this atom? unused ***/
    public double evaluate() {
	
	return 0;
    }


	  public void update() {}

    public String toString() {
	return atom.toString();
    }

    



}
