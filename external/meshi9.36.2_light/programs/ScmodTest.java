package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.solvation.SolvationCreator;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentException;
import meshi.util.*;

import java.io.IOException;

/**
 * Created by chen on 04/05/2015.
 */
public class ScmodTest {
    public static void main(String[] args) throws IOException, UpdateableException, EvaluationException, AlignmentException {
        EnergyCreator[] energyCreators = {
                new BondCreator(),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                new AtomicPairwisePMFSummaCreator(),
                new RamachandranSidechainEnergyCreator() ,
                new SolvationCreator()
        };
        MeshiProgram.initRandom(0);

        AtomList atoms = new AtomList("model.pdb");
        CommandList commands = new CommandList("commands");
        Protein model = new Protein(atoms, ResidueExtendedAtomsCreator.creator);
        model.printAtomsToFile("modelForSPDBV.pdb");

        model = new Protein(atoms, ResidueExtendedAtomsCreator.creator);
        Utils.relax(model, energyCreators, commands);
        model.atoms().molecularSystem.createDistanceMatrix("test");
        Utils.relax(model, energyCreators, commands);
        model.printAtomsToFile("out2.pdb");

        Scmod.scmod(commands,model,10,true, MeshiProgram.randomNumberGenerator());
        Utils.relax(model, energyCreators, commands);
        Scmod.scmod(commands,model,10,true, MeshiProgram.randomNumberGenerator());
        Utils.relax(model, energyCreators, commands);
        model.printAtomsToFile("out1.pdb");


    }
}
