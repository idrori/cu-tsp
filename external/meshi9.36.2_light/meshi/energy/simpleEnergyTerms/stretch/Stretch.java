package meshi.energy.simpleEnergyTerms.stretch;

import meshi.energy.*;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.atoms.Atom;
import meshi.util.UpdateableException;
import meshi.util.Utils;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 30/04/14
 * Time: 23:32
 * To change this template use File | Settings | File Templates.
 */
public class Stretch extends AbstractEnergy{
    private Atom atom1, atom2;
    private double target;

    private EnergyInfoElement info;

    public FreeDistance getDistance() {
        return distance;
    }

    private FreeDistance distance;
    public Stretch(Atom atom1, Atom atom2, double target, EnergyInfoElement energyInfoElement) {
        super(toArray(), energyInfoElement, EnergyType.DIFFERENTIAL);
        System.out.println("Stretch energy "+target+" "+weight+"\n"+atom1+"\n"+atom2);
        this.atom1 = atom1;
        this.atom2 = atom2;
        this.target = target;
        this.info = energyInfoElement;
        distance = new FreeDistance(atom1,atom2);
        comment = "Stretch";
    }

    public void update(int numberOfUpdates) throws UpdateableException {}
    public void test(TotalEnergy energy, Atom atom) {}

    public void evaluateAtoms() {

    }

        public EnergyInfoElement evaluate() {
        distance.update();
//            distance.update();
//            Utils.printDebug(this,"xxxxxxxxxxxxxxxxxxxxxx "+distance.distance);

        //double dis = distance.distance();
        double dis = atom1.x()-atom2.x();
        double diff = dis-target;
        if (diff > 0) {
            info.setValue(0);
            return info;
        }
        else {
            double diff2 = diff*diff;
            double diffPlus1 = diff+1;
            double energy = -diff2*weight/diffPlus1;
            double force = weight*(2*diff/diffPlus1-diff2/(diffPlus1*diffPlus1));
            info.setValue(energy);
            info.getChildren().get(0).setValue(new Double(distance.distance));

            double dx = 1;
//            double dx  = distance.dDistanceDx();
//            double dy  = distance.dDistanceDy();
//            double dz  = distance.dDistanceDz();

            atom1.addToFx(force*dx);
//            atom1.addToFy(force*dy);
//            atom1.addToFz(force*dz);
            atom2.addToFx(-force*dx);
//            atom2.addToFy(-force * dy);
//            atom2.addToFz(-force * dz);
        }

        return info;
    }

    public EnergyElement createElement(Object obj, Parameters parameters) {return null;}

}
