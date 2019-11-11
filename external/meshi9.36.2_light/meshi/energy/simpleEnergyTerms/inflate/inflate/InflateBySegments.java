/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.inflate.inflate;

import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Residue;
import meshi.molecularElements.SegmentList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.MeshiProgram;

import java.util.ArrayList;


public class InflateBySegments extends Inflate {
/*
    private static ArrayList dummy = null;
    private static AtomList atomList, atomListCopy;
    private DistanceMatrix distanceMatrix;
    private Filter filter;
    private double rmsTarget;
    private SegmentList segments;
*/

    public InflateBySegments() {
    }

    public InflateBySegments(DistanceMatrix distanceMatrix, double rmsTarget, EnergyInfoElement info, SegmentList segments) {
        super(distanceMatrix, rmsTarget, info, segments);
        comment = "InflateBySegments ;)";
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
        return true;
    }


    public double weight() {
        return weight;
    }

}