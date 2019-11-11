/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.beta;

import meshi.energy.EnergyElement;
import meshi.energy.beta.peptideBonds.PeptideBond;
import meshi.energy.beta.peptideBonds.PeptideBondDistance;
import meshi.energy.beta.peptideBonds.PeptideBondDistanceList;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.Utils;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 18/01/2010
 * Time: 06:25:43
 */
public class BetaEnergyElement extends EnergyElement {
    public static double ALPHA = 1000;
    public static final double OO_THRESHOLD = 3.2;
    public static final double HH_THRESHOLD = 2.7;
    private PeptideBondDistanceList distances;
    private ArrayList<Diff> diffs1, diffs2;
    private double weight;
    private ArrayList<PeptideBondDistanceList> neighbors;
    private PeptideBondDistanceList neighborDistances;
    private boolean on;
    private PeptideBond myPeptideBond = null;

    public boolean isOn() {
        return on;
    }

    public void on() {
        if ((myPeptideBond != null) & (diffs1 != null))
            on = true;
    }

    public int numberOfNeighbors;

    public BetaEnergyElement(Residue residue, ArrayList<Segment> segments, double weight) {
        if (residue.peptideBond() == null) {
            Utils.println(residue + " has null peptide bond ");
            on = false;
            return;
        }
        myPeptideBond = residue.peptideBond();
        this.weight = weight;
        on = true;
        diffs1 = diffs2 = null;
        Segment mySegment = myPeptideBond.segment();
        distances = new PeptideBondDistanceList(segments.size() * 5);
        atoms = new AtomList(myPeptideBond.H.molecularSystem);
        atoms.add(myPeptideBond.H);
        atoms.add(myPeptideBond.O);
        // Second pass add all the distances of a neighboring segment
        neighbors = myPeptideBond.neighbors();
        numberOfNeighbors = myPeptideBond.numberOfNeighbors();
        for (Segment segment : segments) {
            if (segment != mySegment) {
                neighborDistances = neighbors.get(segment.ID());
                for (Residue res : segment) {
                    if (res.peptideBond() != null) {
                        PeptideBondDistance distance = new PeptideBondDistance(myPeptideBond, res.peptideBond());
                        if (!neighborDistances.isEmpty()) {//That is, there is at least one contact
                            if (distance.distance() > BetaEnergy.TARGET_DISTANCE)
                                neighborDistances.add(distance);// the smaller distances have already entered
                        } else distances.add(distance);
                    }
                }
            }
        }


        if (numberOfNeighbors >= 3) {// This is a bit too much . I will break all these contacts. Hopefully they will rearrange better next time
            System.out.println(this + " has " + numberOfNeighbors + " neighbors.");
            weight = weight * -1;
            distances.clear();
            for (PeptideBondDistanceList pbdl : neighbors)
                for (PeptideBondDistance pbd : pbdl)
                    distances.add(pbd);
            diffs1 = getDiffs(distances);
        } else if (numberOfNeighbors == 2) {// This may be OK if part of a sheet, but very bad as a part of a prism.
            PeptideBondDistanceList pbdl1 = null, pbdl2 = null;
            for (PeptideBondDistanceList pbdl : neighbors) {
                if (!pbdl.isEmpty()) {
                    if (pbdl1 == null) pbdl1 = pbdl;
                    else if (pbdl2 == null) pbdl2 = pbdl;
                    else throw new RuntimeException("This is weird.");
                }
            }
            distances.clear();
            for (PeptideBondDistance pbd1 : pbdl1) {
                for (PeptideBondDistance pbd2 : pbdl2) {
                    if (dis(myPeptideBond.H, pbd1.pb2.O) <= BetaEnergy.TARGET_DISTANCE) {
                        if (dis(myPeptideBond.H, pbd2.pb2.O) <= BetaEnergy.TARGET_DISTANCE) {
                            distances.add(pbd1);
                            distances.add(pbd2);
                        }
                    } else {
                        if (dis(myPeptideBond.O, pbd1.pb2.H) <= BetaEnergy.TARGET_DISTANCE) {
                            if (dis(myPeptideBond.O, pbd2.pb2.H) <= BetaEnergy.TARGET_DISTANCE) {
                                distances.add(pbd1);
                                distances.add(pbd2);
                            }
                        }
                    }
                }
            }
            if (!distances.isEmpty()) {
                System.out.println(this + " is part of a prism. I'll break it.");
                weight = weight * -1;
                diffs1 = getDiffs(distances);
            } else on = false; // Normal situation nothing to build or destroy
        } else if (numberOfNeighbors == 1) {
            for (PeptideBondDistanceList dl : neighbors) {
                if (!dl.isEmpty()) {
                    diffs1 = getDiffs(dl);
                }
            }
            diffs2 = getDiffs(distances, 11);
        } else {
            diffs1 = getDiffs(distances);
        }

        //System.out.println("Created "+this+" "+diffs1.size() );
    }

    public String toString() {
        return "BetaEnergyElement with " + myPeptideBond + " " + on;
    }

    private ArrayList<Diff> getDiffs(PeptideBondDistanceList tempDistances) {
        return getDiffs(tempDistances, 1000);
    }

    private ArrayList<Diff> getDiffs(PeptideBondDistanceList tempDistances, double threshold) {
        ArrayList<Diff> out = new ArrayList<Diff>();
        PeptideBondDistanceList distances = new PeptideBondDistanceList(BetaEnergy.NUMBER_OF_DISTANCES);
        PeptideBondDistanceList sortedList = tempDistances.sort(tempDistances);
        int maxNumberOfDistance = BetaEnergy.NUMBER_OF_DISTANCES;
        if (maxNumberOfDistance > sortedList.size()) maxNumberOfDistance = sortedList.size();
        for (int i = 0; i < maxNumberOfDistance; i++) {
            PeptideBondDistance distance = sortedList.get(i);
            if (distance.distance() < threshold) {
                if (distance.pb2.numberOfNeighbors() < 2) {
                    atoms.add(distance.pb2.H);
                    atoms.add(distance.pb2.O);
                    distances.add(distance);
                }
            }
        }
        double target;
        if (weight > 0) // that is we want to encourage this bond
            target = BetaEnergy.TARGET_DISTANCE;
        else // that is we want to break a contact
            target = 0.1;
        for (PeptideBondDistance d : distances)
            out.add(new Diff(d, target));
        return out;

    }

    public void setAtoms() {
        throw new RuntimeException("Please do not use this method.");
    }


    public double evaluateAtoms() {
        return (evaluateAtoms(diffs1) + evaluateAtoms(diffs2));
    }

    public double evaluateAtoms(ArrayList<Diff> diffs) {
        if (!on) return 0;
        double energy = evaluate();
        if (atoms == null) throw new RuntimeException("getAtoms == null");
        for (Atom atom : atoms)
            atom.addEnergy(energy / atoms.size());
        return energy;
    }

    public double evaluate() {
        double energy = 0;// = clashes.evaluate();
        if (!on) return energy;
        energy += evaluate(diffs1);
        if (diffs2 != null) energy += evaluate(diffs2);
        if (energy == 0) on = false;
        return energy;
    }

    public void off() {
        on = false;
    }

    public double evaluate(ArrayList<Diff> diffs) {
        if (diffs == null) throw new RuntimeException("This is weird4");
        for (Diff diff : diffs) diff.update();
        double tempEnergy = 1, energy, tempEnergyPlusALPHA, tempEnergyPlusALPHA2;
        double dEnergy1DdiffX12, dEnergy1DdiffY12, dEnergy1DdiffZ12, dEnergy1DdiffX21, dEnergy1DdiffY21, dEnergy1DdiffZ21;
        double dEnergyDdiffX12, dEnergyDdiffY12, dEnergyDdiffZ12, dEnergyDdiffX21, dEnergyDdiffY21, dEnergyDdiffZ21;
        for (Diff diff : diffs) {
            tempEnergy *= diff.diff;
        }
        if (tempEnergy == 0) return 0;
        tempEnergyPlusALPHA = tempEnergy + ALPHA;
        tempEnergyPlusALPHA2 = tempEnergyPlusALPHA * tempEnergyPlusALPHA;
        energy = tempEnergy / (tempEnergyPlusALPHA);
        for (Diff diff : diffs) {
            dEnergy1DdiffX12 = diff.dDiffDdX12 * tempEnergy / diff.diff;
            dEnergy1DdiffY12 = diff.dDiffDdY12 * tempEnergy / diff.diff;
            dEnergy1DdiffZ12 = diff.dDiffDdZ12 * tempEnergy / diff.diff;
            dEnergy1DdiffX21 = diff.dDiffDdX21 * tempEnergy / diff.diff;
            dEnergy1DdiffY21 = diff.dDiffDdY21 * tempEnergy / diff.diff;
            dEnergy1DdiffZ21 = diff.dDiffDdZ21 * tempEnergy / diff.diff;
            dEnergyDdiffX12 = dEnergy1DdiffX12 * ALPHA / tempEnergyPlusALPHA2;
            dEnergyDdiffY12 = dEnergy1DdiffY12 * ALPHA / tempEnergyPlusALPHA2;
            dEnergyDdiffZ12 = dEnergy1DdiffZ12 * ALPHA / tempEnergyPlusALPHA2;
            dEnergyDdiffX21 = dEnergy1DdiffX21 * ALPHA / tempEnergyPlusALPHA2;
            dEnergyDdiffY21 = dEnergy1DdiffY21 * ALPHA / tempEnergyPlusALPHA2;
            dEnergyDdiffZ21 = dEnergy1DdiffZ21 * ALPHA / tempEnergyPlusALPHA2;
            diff.atom12H().addToFx(-weight * dEnergyDdiffX12);
            diff.atom12H().addToFy(-weight * dEnergyDdiffY12);
            diff.atom12H().addToFz(-weight * dEnergyDdiffZ12);
            diff.atom12O().addToFx(weight * dEnergyDdiffX12);
            diff.atom12O().addToFy(weight * dEnergyDdiffY12);
            diff.atom12O().addToFz(weight * dEnergyDdiffZ12);
            diff.atom21H().addToFx(-weight * dEnergyDdiffX21);
            diff.atom21H().addToFy(-weight * dEnergyDdiffY21);
            diff.atom21H().addToFz(-weight * dEnergyDdiffZ21);
            diff.atom21O().addToFx(weight * dEnergyDdiffX21);
            diff.atom21O().addToFy(weight * dEnergyDdiffY21);
            diff.atom21O().addToFz(weight * dEnergyDdiffZ21);
        }
        return weight * energy;
    }

    private static class Diff {
        PeptideBondDistance distance;
        double target;
        double d, d2, d3, d4, d3Plus1, d3Plus1Sq, diff;
        double dDiffDd, dDiffDdX12, dDiffDdY12, dDiffDdZ12, dDiffDdX21, dDiffDdY21, dDiffDdZ21;

        public Diff(PeptideBondDistance distance, double target) {
            this.distance = distance;
            this.target = target;
        }

        public void update() {
            distance.update();
            d = distance.distance() - target;
            if (d < 0) {
                diff = 0;
                dDiffDd = 0;
            } else {
                d2 = d * d;
                d3 = d2 * d;
                d4 = d3 * d;
                d3Plus1 = d3 + 1;
                d3Plus1Sq = d3Plus1 * d3Plus1;
                diff = d4 / d3Plus1;
                dDiffDd = 4 * d3 / d3Plus1 - d4 / d3Plus1Sq * 3 * d2;
            }
            dDiffDdX12 = dDiffDd * distance.dDistanceDx12();
            dDiffDdY12 = dDiffDd * distance.dDistanceDy12();
            dDiffDdZ12 = dDiffDd * distance.dDistanceDz12();
            dDiffDdX21 = dDiffDd * distance.dDistanceDx21();
            dDiffDdY21 = dDiffDd * distance.dDistanceDy21();
            dDiffDdZ21 = dDiffDd * distance.dDistanceDz21();
        }

        public Atom atom12H() {
            return distance.pb1.H;
        }

        public Atom atom12O() {
            return distance.pb2.O;
        }

        public Atom atom21H() {
            return distance.pb2.H;
        }

        public Atom atom21O() {
            return distance.pb1.O;
        }
    }


    public double dis(Atom a1, Atom a2) {
        return (new FreeDistance(a1, a2)).distance();
    }
}
