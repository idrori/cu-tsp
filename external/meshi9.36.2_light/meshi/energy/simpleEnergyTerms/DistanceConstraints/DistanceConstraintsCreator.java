package meshi.energy.simpleEnergyTerms.DistanceConstraints;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.SimpleEnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.FreeDistance;
import meshi.geometry.FreeDistanceList;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.Utils;
import meshi.util.info.InfoType;

/**
 *
 */
public class DistanceConstraintsCreator extends SimpleEnergyCreator implements KeyWords {
    public DistanceConstraintsCreator() {
        super(InfoType.DISTANCE_CONSTRAINTS);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,
                                           DistanceMatrix dm, CommandList commands){

        EnergyInfoElement energyInfoElement = new EnergyInfoElement(InfoType.DISTANCE_CONSTRAINTS,
                                                                    "distanced energy", weight);
        CommandList myCommands = null;
        FreeDistanceList distanceList = new FreeDistanceList();
        DistanceConstraintsParametersList targets = new DistanceConstraintsParametersList();
        if (commands.keyExists("constraint"))    {
            myCommands = commands.firstWordFilter("constraint");
            Atom atom1, atom2;
            for (Command command: myCommands){
                atom1 = null;
                atom2 = null;
                int residueNumber1 = command.secondWordInt();
                String atomName1   = command.thirdWord();
                int residueNumber2 = command.fourthWordPositiveInt();
                String atomName2   = command.fifthWord();
                double target   = command.sixthWordDouble();
                Chain chain = protein.chain();
                Residue residue1 = chain.get(residueNumber1);
                for (Atom atom : residue1.getAtoms()) {
                     if (atom.name().equals(atomName1))
                        atom1 = atom;
                }
                if (atom1 == null) throw new RuntimeException("This is weird distanceConstraint "+command);
                Residue residue2 = chain.get(residueNumber2);
                for (Atom atom : residue2.getAtoms()) {
                    if (atom.name().equals(atomName2))
                            atom2 = atom;
                }
                if (atom2 == null)
                    throw new RuntimeException("This is weird distanceConstraint "+command+"\n"+
                                               residue1+" "+atomName1+"    "+residue2+"  "+atomName2);
                distanceList.add(new FreeDistance(atom1,atom2));
                targets.add(new DistanceConstraintsParameters(target));
                Utils.println("Distance constraints   "+target+"\n"+atom1+"\n"+atom2);
            }
        }

        term =  new DistanceConstraints(distanceList, targets,energyInfoElement);
       // System.out.println(term);
        return term;
    }

}
