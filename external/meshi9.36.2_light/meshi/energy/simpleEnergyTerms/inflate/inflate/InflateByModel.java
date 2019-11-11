/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.inflate.inflate;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.Parameters;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.MeshiAttribute;
import meshi.util.MeshiProgram;
import meshi.util.filters.Filter;

import java.io.File;
import java.util.ArrayList;


public class InflateByModel extends Inflate {
//    private static ArrayList dummy = null;
//    private static AtomList atomList, atomListCopy;
//    private DistanceMatrix distanceMatrix;
//    private Filter filter;
//    private double rmsTarget;
//    private SegmentList segments;
    private ArrayList<Segment[]> segmentPairs;
    private File directory;
    private Filter fileFilter;

    public InflateByModel() {
    }

    public InflateByModel(DistanceMatrix distanceMatrix, double rmsTarget, EnergyInfoElement info, SegmentList segments, File directory, Filter fileFilter) {
        super(distanceMatrix, rmsTarget, info, segments);
        comment = "InflatePerSegment ;)";
        this.directory = directory;
        this.fileFilter = fileFilter;
        this.weight /= 10;


    }


    public boolean createElementsList(AtomList atoms) {
        elementsList = new ArrayList();
        Atom atomI, atomJ;
        Residue residueI, residueJ;
        segmentPairs = new ArrayList<Segment[]>();
        Protein otherModel;
        File otherFile;
        Residue otherResidueI, otherResidueJ;
        ArrayList<File> files = new ArrayList<File>();
        for (File file : directory.listFiles()) {
            if (fileFilter.accept(file)) files.add(file);
        }
        if (files.size() == 0) return false;
        otherFile = files.get(MeshiProgram.randomNumberGenerator().nextInt(files.size()));
        otherModel = Protein.getCAproteinFromApdbFile(otherFile);
        int counter = 0;
        FreeDistance distance, otherDistance;

        for (int i = 0; i < atoms.size(); i++) {
            atomI = atoms.get(i);
            if (atomI.backboneCA()) {
                residueI = atomI.residue();
                otherResidueI = getOtherResidue(otherModel, residueI);
                if (otherResidueI != null) {
                    for (int j = 0; j < i; j++) {
                        atomJ = atoms.get(j);
                        if (atomJ.backboneCA()) {
                            residueJ = atomJ.residue();
                            otherResidueJ = getOtherResidue(otherModel, residueJ);
                            if (otherResidueJ != null) {
                                distance = new FreeDistance(atomI, atomJ);
                                otherDistance = new FreeDistance(otherResidueI.ca(), otherResidueJ.ca());
                                if (Math.abs(distance.distance() - otherDistance.distance()) > 1) {
                                    elementsList().add(new InflateEnergyElement(distance, otherDistance.distance, 0, weight));
                                    counter++;
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.println("Extracted " + counter + " distance constraints from " + otherFile);
        return true;
    }


    private static Residue getOtherResidue(Protein model, Residue myResidue) {
        Residue out = null;
        int number = myResidue.number();
        if (number < model.chain().size()) {
            out = model.chain().get(number);
            if ((!out.dummy()) &&
                    (out.type() == myResidue.type())) return out;
        }
        return null;
    }

    private boolean uniqueSegmentPair(Residue residueI, Residue residueJ) {
        Segment segmentI = ((SegmentAttribute) residueI.getAttribute(MeshiAttribute.SEGMENT_ATTRIBUTE)).segment();
        Segment segmentJ = ((SegmentAttribute) residueJ.getAttribute(MeshiAttribute.SEGMENT_ATTRIBUTE)).segment();
        if (segmentI == segmentJ) return false;
        for (Segment[] segmentPair : segmentPairs) {
            if ((segmentPair[0] == segmentI) & (segmentPair[1] == segmentJ)) return false;
            if ((segmentPair[1] == segmentI) & (segmentPair[0] == segmentJ)) return false;
        }
        Segment[] newPair = new Segment[2];
        newPair[0] = segmentI;
        newPair[1] = segmentJ;
        segmentPairs.add(newPair);
        return true;
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        throw new RuntimeException("Do not use this");
    }


    public double weight() {
        return weight;
    }

}