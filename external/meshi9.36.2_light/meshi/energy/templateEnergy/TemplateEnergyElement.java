package meshi.energy.templateEnergy;
import meshi.energy.EnergyElement;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;


public class TemplateEnergyElement extends EnergyElement {
    private FreeDistance distance;
    private Atom movingAtom;
    private double weight;
    private double EPSILON = 0.01;
    private double BETA = 1000;
    private boolean on = true;
    private double dis, disPlusEpsilon,disDivaidedByDisPlusEpsilon;
    private double g, dGdD, energy, dEdD;
    private double gDividedByGplusBETA;

    public TemplateEnergyElement(Atom frozenAtom, Atom movingAtom, double weight) {
	distance = new FreeDistance(movingAtom, frozenAtom);
	this.movingAtom = movingAtom;
	setAtoms();
	this.weight = weight;
    }
    public void update() {if (on) distance.update();}

    public double evaluate() {
	if (!on) return 0;
	dis = distance.distance();
	if (dis == 0) return 0; 
	disPlusEpsilon = dis+EPSILON;
	disDivaidedByDisPlusEpsilon = dis/disPlusEpsilon;
	g = weight*dis*disDivaidedByDisPlusEpsilon;
	dGdD = weight*disDivaidedByDisPlusEpsilon*(2-disDivaidedByDisPlusEpsilon);
	gDividedByGplusBETA = g/(g+BETA);
        energy = BETA*gDividedByGplusBETA; 
	dEdD = BETA*dGdD*(1-gDividedByGplusBETA)/(g+BETA);
	movingAtom.addToFx(-dEdD * distance.dDistanceDx());	
	movingAtom.addToFy(-dEdD * distance.dDistanceDy());	
	movingAtom.addToFz(-dEdD * distance.dDistanceDz());
	return energy;
    }

    protected void setAtoms(){
	atoms = new AtomList(movingAtom.molecularSystem);
	atoms.add(movingAtom);
    }
    public Residue residue() {
	return movingAtom.residue();
    }
    public void off() {on = false;}
    public void offIfUnhappy(double threshold) {
	if (on) {
	    update();
	    if (distance.distance() >= threshold) on = false;
	}
    }
}
