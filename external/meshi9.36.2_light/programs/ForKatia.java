package programs;

import meshi.energy.EvaluationException;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.Loop;
import meshi.sequences.AlignmentException;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.UpdateableException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by chen on 04/05/2015.
 */
public class ForKatia  {
    public static void main(String[] args) throws IOException, UpdateableException, EvaluationException, AlignmentException{
        MeshiProgram.initRandom(0);
        AtomList atoms = new AtomList("model.pdb");
        CommandList commands = new CommandList("commands");
        Protein model = new Protein(atoms, ResidueExtendedAtomsCreator.creator);
        model.printAtomsToFile("modelForSPDBV.pdb");

        List<Loop> loops = new ArrayList();
        Loop.getLoops(model, loops);
        Loop.addLoopsQuickAndDirty(model, loops, commands, 4, 1000,2);

        model.printAtomsToFile("out.pdb");
    }
}
