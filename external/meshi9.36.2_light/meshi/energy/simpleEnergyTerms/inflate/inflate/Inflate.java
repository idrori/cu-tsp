/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.inflate.inflate;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Residue;
import meshi.molecularElements.Segment;
import meshi.molecularElements.SegmentAttribute;
import meshi.molecularElements.SegmentList;
import meshi.molecularElements.atoms.*;
import meshi.optimizers.Minimizer;
import meshi.util.MeshiAttribute;
import meshi.util.MeshiProgram;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;

import java.util.ArrayList;


public class Inflate extends SimpleEnergyTerm {
    private static int debugCounter = 0;
//    private static ArrayList dummy = null;
    private AtomList atomList, atomListCopy;
    private DistanceMatrix distanceMatrix;
//    private Filter filter;
    private double rmsTarget;
    private SegmentList segments;

    public Inflate() {
    }

    public Inflate(DistanceMatrix distanceMatrix, double rmsTarget, EnergyInfoElement info, SegmentList segments) {
        super(toArray(distanceMatrix), null, info);
        comment = "Inflate ;)";
        this.distanceMatrix = distanceMatrix;
        off();
        this.rmsTarget = rmsTarget;
        this.segments = segments;

    }

    public EnergyInfoElement evaluate() throws EvaluationException{
        boolean targetReached;
        try {
            targetReached = atomList.getRms(atomListCopy) > rmsTarget;
        }
        catch (Exception ex) {
            Utils.println("inflate failed due to:\n" + ex + "\n" + "Trying to continue.");
            targetReached = true;
        }
        if (targetReached) {
            Utils.println("inflate reached RMS of " + rmsTarget);
            Minimizer.terminator.kill("InflateBySegments reached RMS of " + rmsTarget);
        }
        return super.evaluate();
    }

    public void off() {
        super.off();
    }

   public void reset() {
       MolecularSystem copyMolecularSystem = new MolecularSystem();
        Utils.println(this+"resetting");
        atomList = new AtomList(distanceMatrix.molecularSystem);
        for (AtomCore atom : distanceMatrix.molecularSystem)
            atomList.add(atom.atom);
        atomListCopy = Utils.duplicateInAnewMolecularSystem(atomList,copyMolecularSystem);
        createElementsList(atomList);
    }


    public boolean areWeThere() {
        boolean targetReached;
        double rms = atomList.getRms(atomListCopy);
        Utils.println("Inflate rms "+rms + " target "+rmsTarget);
        try {
            targetReached = rms > rmsTarget;
        }
        catch (Exception ex) {
            System.out.println("inflate failed due to:\n" + ex + "\n" + "Trying to continue.");
            targetReached = false;
        }
        if (targetReached) {
            System.out.println("inflate reached RMS of " + rmsTarget);


            super.off();
            return true;
        }
        return false;
    }


    public boolean restart() {
        Utils.println(this+"restarting weight = "+weight);
//        if (!on) return false;
         boolean areWeThere = areWeThere();
         if (!areWeThere)   {
             Utils.println("Restarting "+this+"     Rms is "+atomList.getRms(atomListCopy)+" target is"+rmsTarget);
            createElementsList(atomList);
         }
        return (! areWeThere);
    }

    public void scaleWeight(double factor){
        super.scaleWeight(factor);
        for (Object element : elementsList())
            ((InflateEnergyElement)element).scaleWeight(factor);
    }

    public boolean createElementsList(AtomList atoms) {
        elementsList = new ArrayList();
        Atom atomI, atomJ;
        Residue residueI, residueJ;

        for (int i = 0; i < atoms.size(); i++) {
            atomI = atoms.get(i);
            if (atomI.backboneCA()) {
                residueI = atomI.residue();
                if (MeshiProgram.randomNumberGenerator().nextDouble() < 0.3) {
                    for (int j = 0; j < i; j++) {
                        atomJ = atoms.get(j);
                        if (atomJ.backboneCA()) {
                            residueJ = atomJ.residue();
                            if (notOnTheSameSSfragment(residueI, residueJ)) {
                                if ((new FreeDistance(atomI, atomJ)).distance() < 8) {
                                    if (MeshiProgram.randomNumberGenerator().nextDouble() < 0.3) {
                                        elementsList().add(new InflateEnergyElement(new FreeDistance(atomI, atomJ), weight));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return true;
    }


    private boolean notOnTheSameSSfragment(Residue residueI, Residue residueJ) {
        Segment segmentI = ((SegmentAttribute) residueI.getAttribute(MeshiAttribute.SEGMENT_ATTRIBUTE)).segment();
        Segment segmentJ = ((SegmentAttribute) residueJ.getAttribute(MeshiAttribute.SEGMENT_ATTRIBUTE)).segment();
        return (segmentI != segmentJ);
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        throw new RuntimeException("Do not use this");
    }

    public static class TargetFilter implements Filter {
        private double target;

        public TargetFilter(double target) {
            this.target = target;
        }

        public boolean accept(Object obj) {
            AtomPair ap = (AtomPair) obj;
            if (ap.atom1().distanceFrom(ap.atom2()) < target) return true;
            return false;
        }
    }

    public double weight() {
        return weight;
    }

}