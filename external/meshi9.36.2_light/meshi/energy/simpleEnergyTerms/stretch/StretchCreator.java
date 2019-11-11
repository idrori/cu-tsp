package meshi.energy.simpleEnergyTerms.stretch;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 30/04/14
 * Time: 23:21
 * To change this template use File | Settings | File Templates.
 */
public class StretchCreator extends EnergyCreator{
    private double target;
    private Atom atom1, atom2;
    public StretchCreator() {
        super(InfoType.STRETCH);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix dm, CommandList commands){

        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.STRETCH, "stretch", weight);
        energyInfoElement.getChildren().add(new EnergyInfoElement(InfoType.STRETCH,"endToEnd"));
        atom1 = null;
        for (Atom atom: protein.atoms()){
            if (!atom.nowhere()) {
                if (atom1 == null) atom1 = atom;
                atom2 = atom;
            }
        }
        target = 3.8*(atom2.residueNumber()-atom1.residueNumber());
        term =  new Stretch(atom1,atom2,target,energyInfoElement);
        System.out.println(term);
        return term;

        }



    }
