package meshi.energy.simpleEnergyTerms.stretch;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomPair;
import meshi.util.UpdateableException;

import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 30/04/14
 * Time: 23:32
 * To change this template use File | Settings | File Templates.
 */
public class StretchAll extends SimpleEnergyTerm{

    private EnergyInfoElement info;
    private int modelLength;


    public StretchAll(Protein model, EnergyInfoElement energyInfoElement) {
        super(toArray(), null, energyInfoElement);
        this.info = energyInfoElement;
        comment = "StretchAll";
        ArrayList<AtomPair> atomPairs = new ArrayList();
        Atom[] CAs = model.chain().atoms().CAFilter().toArrayOfAtoms();
        for (int iAtom = 0; iAtom < CAs.length/2 - 20; iAtom++) {
            int jAtom = CAs.length-1-iAtom;
            atomPairs.add(new AtomPair(CAs[iAtom], CAs[jAtom]));
        }
        modelLength = model.chain().size();
        createElementsList(atomPairs);
    }

    public void update(int numberOfUpdates) throws UpdateableException {}
    public void test(TotalEnergy energy, Atom atom) {}

    public void evaluateAtoms() {

    }


    public void createElementsList(ArrayList baseList) {
        elementsList = new ArrayList();

        for (Object baseElement : baseList) {
            EnergyElement newElement = createElement(baseElement,null);
            elementsList.add(newElement);
        }
    }

    public EnergyElement createElement(Object obj, Parameters parameters) {
        return new StretchAllElement((AtomPair) obj,-10.0/modelLength);
    }

    private static class StretchAllElement extends EnergyElement {
        private Atom atom1, atom2;
        private double startDistance;
        private FreeDistance distance;
        private double weight;

        public StretchAllElement(AtomPair pair, double weight) {
            atom1 = pair.atom1();
            atom2 = pair.atom2();
            distance = new FreeDistance(atom1,atom2);
            startDistance = distance.distance()*2;
            this.weight = weight;
        }

        public void setAtoms() {
            atoms.add(atom1);
            atoms.add(atom2);
        }

        public double evaluate() {
            distance.update();

            double dis = distance.distance();
            if (dis < startDistance) return 0;
            double diff = dis-startDistance;
            double diff2 = diff*diff;
            double diffPlusOne = diff + 1;
            double dx  = distance.dDistanceDx();
            double dy  = distance.dDistanceDy();
            double dz  = distance.dDistanceDz();
            double energy = weight*diff2/diffPlusOne;
            double dEnergyDdiff = weight*(2*diff/diffPlusOne-diff2/(diffPlusOne*diffPlusOne));
            atom1.addToFx(-dEnergyDdiff * dx);
            atom1.addToFy(-dEnergyDdiff * dy);
            atom1.addToFz(-dEnergyDdiff * dz);
            atom2.addToFx(dEnergyDdiff * dx);
            atom2.addToFy(dEnergyDdiff * dy);
            atom2.addToFz(dEnergyDdiff*dz);
            return energy;
        }
    }

}
