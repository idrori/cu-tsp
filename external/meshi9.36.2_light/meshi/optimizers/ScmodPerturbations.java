/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.EvaluationException;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ResidueList;
import meshi.util.CommandList;
import meshi.util.Scmod;
import meshi.util.UpdateableException;
import meshi.util.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 24/02/2010
 * Time: 23:09:23
 * To change this template use File | Settings | File Templates.
 */
public class ScmodPerturbations implements Perturbation {
    private Protein model;
    private CommandList commands;
    private ResidueList conservedResidues;

    public ScmodPerturbations(Protein model, CommandList commands, ResidueList conservedResidues) {
        this.model = model;
        this.commands = commands;
        this.conservedResidues = conservedResidues;
    }

    public void perturb() throws UpdateableException, EvaluationException{
        Utils.println("\nPerturbation by SC rearrangement\n");
        conservedResidues.freeze();
        Scmod.scmod(commands, model, 4, 50);
     }
    public void restart() {
        throw new RuntimeException("Not implemented yet");
    }
    public void reset() {
        throw new RuntimeException("Not implemented yet");
    }
}
