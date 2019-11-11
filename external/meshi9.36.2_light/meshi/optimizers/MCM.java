/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.AbstractEnergy;
import meshi.energy.EvaluationException;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeSumma5PolarTypesEnergy;
import meshi.sequences.AlignmentException;
import meshi.util.UpdateableException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVol;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.scoringFunctions.Score;
import meshi.util.Logger;
import meshi.util.MeshiProgram;
import meshi.util.Utils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class MCM extends Optimizer {
    protected TotalEnergy energy;
     Minimizer minimizer;
    private TemperatureGenerator temperatureGenerator;
    protected Perturbation perturbation;
    protected Logger log = null;
    protected TetherEnergy allAtomTether = null;
    protected ArrayList<Score> scoreFunctions;
    protected Score optimizationScore;

    public  enum mcmStepResult {
        FAILED, SUCCEEDED, AFTER_PERTURBATION;
        private  int lastSuccess;
        protected static void setLastSuccess(int lastSuccess) {
            FAILED.lastSuccess = lastSuccess;
            SUCCEEDED.lastSuccess = lastSuccess;
            AFTER_PERTURBATION.lastSuccess = lastSuccess;
        }
        public int lastSuccess() {
            return lastSuccess;
        }
    }

    public static enum mcmMode {RELAXATION, OPTIMIZATION}

    public MCM(TotalEnergy energy, ArrayList<Score> scoreFunctions, Score optimizationScore, Minimizer minimizer, Perturbation perturbation,
               TemperatureGenerator temperatureGenerator, int maxSteps) {
        super(energy, maxSteps, 1);
        this.energy = energy;
        this.minimizer = minimizer;
        this.temperatureGenerator = temperatureGenerator;
        this.perturbation = perturbation;
        this.scoreFunctions = scoreFunctions;
        this.optimizationScore = optimizationScore;
        for (AbstractEnergy e : energy.energyTerms()) {
            if (e instanceof TetherEnergy) {
                if (((TetherEnergy) e).allAtomsTether()) allAtomTether = (TetherEnergy) e;
            }
        }
        if (allAtomTether == null)
            throw new RuntimeException("Please add an all getAtoms Tether term for numerical stability at the end of the perturbation");
    }


    public OptimizerStatus run(Logger log)  {
        this.log = log;
        try {
            return run();
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    public OptimizerStatus run() throws OptimizerException, UpdateableException, EvaluationException, AlignmentException, IOException {
        double oldScore;
        double[][] oldCoordinates;
        double currentScore;
        double dS;

        Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();

        energy.setCurrent();
        energy.on();
        allAtomTether.off();

        oldScore = -99999999; // first step always accepted
        Utils.println(" Old score" + oldScore);
        for (int step = 1; step <= maxSteps; step++) {
            energy.setCurrent();
            energy.on();
            ExcludedVol ev = (ExcludedVol)energy.getEnergyTerm(new ExcludedVol());
            ev.off();
            allAtomTether.off();
            double temperature = temperatureGenerator.next();
            energy.evaluate();
            Utils.println("\n step # " + step + "; temperature " + temperature + " \ncurrent energy = \n" + energy.report(step));

            oldCoordinates = getOldCoordinates(energy);
            try {
                if (step > 1) {
                    perturb(perturbation, step);
                    energy.getDistanceMatrix().terminator.reset();
                   // energy.getDistanceMatrix().reset();
                }
                energy.setCurrent();
                energy.on();
                log.mcm(scoreFunctions, energy, step, mcmStepResult.AFTER_PERTURBATION);
                allAtomTether.reset();
                energy.getEnergyTerm(new CooperativeSumma5PolarTypesEnergy()).off();
                minimizer.run();
                allAtomTether.off();
                energy.getEnergyTerm(new CooperativeSumma5PolarTypesEnergy()).on();
                energy.evaluate();
                OptimizerStatus os = minimizer.run();

                Utils.println("MCM step " + os);
                energy.on();
                allAtomTether.off();
                currentScore = ((Double) optimizationScore.score(energy.energyInfo()).getValue()).doubleValue();
                Utils.println("oldScore = " + oldScore + "\n" + "currentScore = " + currentScore + "\n" + energy.report(-999999));
                dS = currentScore - oldScore;
                if (dS < 0) {
                    double rnd = randomNumberGenerator.nextDouble();
                    if (rnd > Math.exp(dS / temperature)) {//That is Metropolis criterion failed
                         if (log != null) {
                                log.mcm(scoreFunctions, energy, step, mcmStepResult.FAILED);
                        }
                        setCoordinates(energy, oldCoordinates);
                        Utils.println("This step failed: dS = " + dS +
                                "  dS/temperature = " + dS / temperature +
                                "  Math.exp(dS/temperature) = " + Math.exp(dS / temperature) +
                                "  rnd = " + rnd);

                    } else {
                         if (log != null) {
                                mcmStepResult.setLastSuccess(step);
                                log.mcm(scoreFunctions, energy, step, mcmStepResult.SUCCEEDED);
                                //energy.distanceMatrix().writeNonBondedList(step+"_SUCCEEDED_1");
                        }
                        Utils.println("This step ded: dS = " + dS +
                                "  dS/temperature = " + dS / temperature +
                                "  Math.exp(dS/temperature) = " + Math.exp(dS / temperature) +
                                "  rnd = " + rnd);
                        oldScore = currentScore;
                    }
                } else {
                    Utils.println("This step succeeded: dS = " + dS);
                     if (log != null) {
                                mcmStepResult.setLastSuccess(step);
                                //  energy.distanceMatrix().writeNonBondedList(step+"_SUCCEEDED_2");
                                log.mcm(scoreFunctions, energy, step, mcmStepResult.SUCCEEDED);
                     }
                    oldScore = currentScore;
                }
            }
            catch (OptimizerException ex) {
                if (Utils.verbose()) energy.test();
                System.out.println("This step failed due to "+ex);
                setCoordinates(energy, oldCoordinates);
                throw new RuntimeException(ex);            }
        }
        energy.on();

        return OptimizerStatus.DONE;
    }

    public static double[][] getOldCoordinates(TotalEnergy energy) {
        double[][] coordinates = energy.coordinates();
        int length = coordinates.length;
        double[][] out = new double[length][2];

        for (int i = 0; i < length; i++) {
            out[i][0] = coordinates[i][0];
            out[i][1] = coordinates[i][1];
        }

        return out;
    }

    public static void setCoordinates(TotalEnergy energy, double[][] toSet) {
        double[][] coordinates = energy.coordinates();
        int length = coordinates.length;
        if (length != toSet.length) throw new RuntimeException("Weird parameters to MCM.setCoordinates");

        for (int i = 0; i < length; i++) {
            coordinates[i][0] = toSet[i][0];
            coordinates[i][1] = toSet[i][1];
        }
    }

    private static Perturbation[] toArray(Perturbation p) {
        Perturbation[] out = {p};
        return out;
    }

    protected void perturb(Perturbation perturbation, int step) throws OptimizerException, EvaluationException, UpdateableException {
         Minimizer.terminator.reset();
         perturbation.reset();
         AtomList atomListCopy = Utils.duplicateInAnewMolecularSystem(energy.atomList(), new MolecularSystem());
                    int counter = 0;
                    try {
                        perturbation.perturb();
                    }
                    catch (OptimizerException ox) {
                        System.out.println("A problem in perturbation #  " + counter);
                        System.out.println(ox);
                        ox.printStackTrace();
                        System.out.println("RMS reached\t" + counter + "\t" + step + "\t" + energy.atomList().getRms(atomListCopy));
                        if (ox.energy() != null) ox.energy().test();
                        System.out.println("\nContinuing\n");
                    }
        Utils.println("RMS perturbation\t" + step + "\t" + energy.atomList().getRms(atomListCopy) + "\t" + log.rms());
        Minimizer.terminator.reset();
       }

}

/*
        InflateBySegments inflate = (InflateBySegments) energy.getEnergyTerm(new InflateBySegments());
        if (inflate == null) throw new RuntimeException("No point in MCM without inflate");
        if (iteration != 1) {
        inflate.isOn();
        System.out.println(optimizer.run());
        }
        inflate.off();
        System.out.println(optimizer.run());
*/
