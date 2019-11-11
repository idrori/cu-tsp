package meshi.optimizers;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVol;
import meshi.molecularElements.atoms.AtomList;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.util.UpdateableException;
import meshi.util.Utils;

import java.io.IOException;
import java.util.ArrayList;

/**
 *
 */
public class Relaxation extends MCM {
    public Relaxation(TotalEnergy energy, ArrayList<Score> scoreFunctions, Score optimizationScore, Minimizer minimizer, Perturbation perturbation,
                      TemperatureGenerator temperatureGenerator, int maxSteps) {
            super(energy, scoreFunctions, optimizationScore, minimizer, perturbation, temperatureGenerator, maxSteps);
    }
    public OptimizerStatus run() throws OptimizerException, UpdateableException,EvaluationException,IOException, AlignmentException {
            double oldScore;
            double[][] oldCoordinates;
            double currentScore;
            double dS;

            oldScore = -9999999; // first step always accepted
            Utils.println(" Old score" + oldScore);
            for (int step = 1; step <= maxSteps; step++) {
                energy.setCurrent();
                energy.on();
                allAtomTether.off();
                energy.evaluate();
                Utils.println("\n relaxation step # " + step + ";   \ncurrent energy = \n" + energy.report(step));

                oldCoordinates = getOldCoordinates(energy);
                try {
                    if (step > 1) perturb(perturbation,step);
                    try {
                        energy.setCurrent();
                    } catch (Exception ex) {
                        System.out.println(" This step failed.\ncould not setResidue the minimization energy to be the current one due to " + ex);
                        setCoordinates(energy, oldCoordinates);
                        continue;
                    }
                    energy.on();
                    ExcludedVol ev = (ExcludedVol)energy.getEnergyTerm(new ExcludedVol());
                    ev.off();

                    // The great moment
                    allAtomTether.reset();
                    energy.evaluate();
                    try {
                        minimizer.run();
                    } catch (Exception ex) {
                        System.out.println(" This step failed.\nFirst (tethered) minimization failed due to " + ex);
                        ex.printStackTrace();
                        setCoordinates(energy, oldCoordinates);
                        continue;
                    }
                    allAtomTether.off();
                    OptimizerStatus os = minimizer.run();

                    System.out.println("MCM step " + os);

                    energy.on();
                    allAtomTether.off();

                    currentScore = ((Double) optimizationScore.score(energy.energyInfo()).getValue()).doubleValue();
                    Utils.println("oldScore = " + oldScore + "\n" + "currentScore = " + currentScore + "\n" + energy.report(999999));
                    dS = currentScore - oldScore;
                    if (dS < 0) {
                             if (log != null) {
                                    log.mcm(scoreFunctions, energy, step, mcmStepResult.FAILED);


                            }
                            setCoordinates(energy, oldCoordinates);
                            Utils.println("This step failed: dS = " + dS +
                                    "  dS = " + dS);

                        } else {
                             if (log != null) {
                                    mcmStepResult.setLastSuccess(step);
                                    log.mcm(scoreFunctions, energy, step, mcmStepResult.SUCCEEDED);
                            }
                            Utils.println("This step succeeded: dS = " + dS +
                                    "  dS = " + dS);
                            oldScore = currentScore;
                                    if (step > 1) return OptimizerStatus.DONE;
                        }
                }
                catch (OptimizerException ox) {
                    energy.test();
                    System.out.println("This step failed due to OptimizerException");
                    setCoordinates(energy, oldCoordinates);
                }


            }
            energy.on();

            return OptimizerStatus.CONVERGED;
        }

}
