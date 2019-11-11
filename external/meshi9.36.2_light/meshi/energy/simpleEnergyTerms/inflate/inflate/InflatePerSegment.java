/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.inflate.inflate;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.Parameters;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Residue;
import meshi.molecularElements.Segment;
import meshi.molecularElements.SegmentAttribute;
import meshi.molecularElements.SegmentList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.MeshiAttribute;
import meshi.util.MeshiProgram;

import java.util.ArrayList;


public class InflatePerSegment extends Inflate {
/*
    private static ArrayList dummy = null;
    private static AtomList atomList, atomListCopy;
    private DistanceMatrix distanceMatrix;
    private Filter filter;
    private double rmsTarget;
    private SegmentList segments;
*/
    private ArrayList<Segment[]> segmentPairs;

    public InflatePerSegment() {
    }

    public InflatePerSegment(DistanceMatrix distanceMatrix, double rmsTarget, EnergyInfoElement info, SegmentList segments) {
        super(distanceMatrix, rmsTarget, info, segments);
        comment = "InflatePerSegment ;)";
        this.weight *= 10;
    }


    public boolean createElementsList(AtomList atoms) {
        elementsList = new ArrayList();
        Atom atomI, atomJ;
        Residue residueI, residueJ;
        segmentPairs = new ArrayList<Segment[]>();

        for (int i = 0; i < atoms.size(); i++) {
            atomI = atoms.get(i);
            if (atomI.backboneCA()) {
                residueI = atomI.residue();
                if (MeshiProgram.randomNumberGenerator().nextDouble() < 0.3) {
                    for (int j = 0; j < i; j++) {
                        atomJ = atoms.get(j);
                        if (atomJ.backboneCA()) {
                            residueJ = atomJ.residue();
                            if (uniqueSegmentPair(residueI, residueJ)) {
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