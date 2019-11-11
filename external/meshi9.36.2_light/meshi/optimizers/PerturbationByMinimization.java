/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.beta.BetaEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleEnergyTerms.inflate.inflate.Inflate;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateByModel;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateBySegments;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflatePerSegment;
import meshi.energy.twoTorsions.TwoTorsionsEnergy;
import meshi.util.*;
import meshi.molecularElements.*;
import meshi.energy.*;

import java.util.*;


public class PerturbationByMinimization implements KeyWords, Perturbation {
    private TotalEnergy energy;
    private Inflate perturbationTerm;
    private Minimizer minimizer;
    private CommandList commands;
    private Protein protein;
    private Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();
    ResidueList conservedResidues;
    private String comment;


    public PerturbationByMinimization(TotalEnergy energy, CommandList commands, Protein protein,
                                                                                      ResidueList conservedResidues,
                                                                                      String comment) throws UpdateableException {
        this.energy = energy;
        perturbationTerm = (InflateBySegments) energy.getEnergyTerm(new InflateBySegments());
        if (perturbationTerm == null)
            perturbationTerm = (InflatePerSegment) energy.getEnergyTerm(new InflatePerSegment());
        if (perturbationTerm == null) perturbationTerm = (Inflate) energy.getEnergyTerm(new Inflate());
        if (perturbationTerm == null) perturbationTerm = (Inflate) energy.getEnergyTerm(new InflateByModel());
        if (perturbationTerm != null) Utils.println("Perturbation term " + perturbationTerm.comment());
        else throw new RuntimeException("Failed to find perturbation term");
        this.commands = commands;
        this.protein = protein;
        minimizer = Utils.getLBFGS(energy, commands, MCM_PERTURBATION);
        this.conservedResidues = conservedResidues;
        this.comment = comment;
    }


    public void reset(){
        Utils.println(this+" resetting");
        perturbationTerm.reset();
    }
    public void restart(){
        Utils.println(this+" restarting");
        perturbationTerm.restart();
    }

    public void perturb() throws OptimizerException, UpdateableException, EvaluationException {
        Utils.println("\nPerturbation by minimization\n");

        energy.setCurrent();
        BetaEnergy beta = (BetaEnergy) energy.getEnergyTerm(new BetaEnergy());
        if (beta != null) beta.reset();
        RamachandranEnergy.resetRamachandran(energy);
        AbstractEnergy flatRamach = energy.getEnergyTerm(new TwoTorsionsEnergy());
        flatRamach.scaleWeight(100);
        energy.evaluate();
        Utils.println("\n Initial perturbation Energy  energy = \n" + energy.report(0));
                for (int j = 0; (j< 3) &(!perturbationTerm.areWeThere()); j++) {
                        perturbationTerm.scaleWeight((j+1)/3.0);
                        Utils.print("PerturbationByMinimization # :"+j+" with "+perturbationTerm+" "+perturbationTerm.weight());
                        minimizer.run(false);
                        if (!perturbationTerm.areWeThere())
                            perturbationTerm.scaleWeight(3.0/(j+1));
                    }
        flatRamach.scaleWeight(0.01);
        Scmod.scmod(commands, protein, 20, true, randomNumberGenerator);
        energy.distanceMatrix().terminator.reset();
    }

    public String toString() {
        return "PerturbationByMinimization " + comment;
    }
}



