/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.EvaluationException;
import meshi.util.UpdateableException;
import meshi.util.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 22/09/2010
 * Time: 09:18:22
 * To change this template use File | Settings | File Templates.
 */
public class CombinedPerturbation implements Perturbation {
    private final Perturbation[] perturbations;
    public static int counter, numberOfPerturbations;

    public CombinedPerturbation(Perturbation[] perturbations) {
        this.perturbations = perturbations;
        counter = 0;
        numberOfPerturbations = perturbations.length;
    }

    public void reset(){
        for (Perturbation perturbation : perturbations){
            Utils.println(this+" resetting");
            perturbation.reset();
        }
    }
    public void restart(){
           Utils.println(this+" restarting");
           for (Perturbation perturbation : perturbations){
               perturbation.restart();
           }
       }

    public void perturb() throws OptimizerException, UpdateableException, EvaluationException {
        for ( int i = 0; (i < 10) & (!Minimizer.terminator.dead()); i++) {
            int perturbationNumber = i % numberOfPerturbations;
            Perturbation perturbation = perturbations[perturbationNumber];
            System.out.println("Combined perturbation " + " " + i + " performed by element number " + perturbationNumber + " " + perturbation);
            perturbation.restart();
            perturbation.perturb();
        }
    }
}
