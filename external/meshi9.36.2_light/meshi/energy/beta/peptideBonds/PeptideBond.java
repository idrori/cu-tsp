/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.beta.peptideBonds;


import meshi.energy.beta.BetaEnergy;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;
import meshi.parameters.*;
import meshi.util.MeshiAttribute;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 22/01/2010
 * Time: 16:22:14
 * To change this template use File | Settings | File Templates.
 */
public class PeptideBond {
    public final Atom H, O;
    public final Residue residue;
    ArrayList<PeptideBondDistanceList> neighbors = null;

    public PeptideBond(Atom H, Atom O) {
        this.H = H;
        this.O = O;
        this.residue = H.residue();
        if (O.residue().number() != residue.number() - 1)
            throw new RuntimeException("Peptide bond is between successive residues.\n" +
                    " Here the residues are " + H.residue() + " and " + O.residue());
        if (!H.type().backboneH()) throw new RuntimeException("Weird peptide bond H atom is " + H);
        if (!O.type().backboneO()) throw new RuntimeException("Weird peptide bond O atom is " + O);
        if (residue.type() == ResidueType.PRO) {
            throw new RuntimeException("Currently peptide bonds of Proline are not supported ");
        }
        if (residue.peptideBond() == null) residue.setPeptideBond(this);
        else if (residue.peptideBond() != this)
            throw new RuntimeException("Only one PeptideBond instance may be associated with a residue");
    }

    public String toString() {
        return "PeptideBond of " + residue + "  " + residue.getAttribute(MeshiAttribute.SEGMENT_ATTRIBUTE);
    }

    public static void attachPeptideBonds(Protein protein) {
        for (Chain chain : protein.chains()) {
            attachPeptideBonds(chain);
        }
    }

    public static void attachPeptideBonds(Chain chain) {
        for (Residue residue : chain) {
            if (residue.dummy()) continue;
            if (residue.type() == ResidueType.PRO) continue;
            if (residue.peptideBond() != null) continue;
            Residue prevResidue = chain.get(residue.number() - 1);
            if (prevResidue.dummy()) continue;
            Atom O = prevResidue.carbonylO();
            Atom H = residue.amideH();
            if (O == null) throw new RuntimeException("weird prev residue for PeptideBond " + prevResidue);
            if (H == null) throw new RuntimeException("Weird residue for PeptideBond " + residue);
            new PeptideBond(H, O); // The new object is attached to the residue in the constructor.
        }
    }

    public Segment segment() {
        return ((SegmentAttribute) residue.getAttribute(MeshiAttribute.SEGMENT_ATTRIBUTE)).segment();
    }


    public void neighbors(ArrayList<Segment> segments) {
        neighbors = new ArrayList<PeptideBondDistanceList>();

        for (Segment seg : segments)
            neighbors.add(new PeptideBondDistanceList(10));

        //First pass - find neighbors
        Segment mySegment = segment();
        PeptideBondDistanceList neighborDistances;
        for (Segment segment : segments) {
            if (segment != mySegment) {
                neighborDistances = neighbors.get(segment.ID());
                for (Residue res : segment) { //
                    if (res.peptideBond() != null) {
                        PeptideBondDistance distance = new PeptideBondDistance(this, res.peptideBond());
                        if (distance.distance() <= BetaEnergy.TARGET_DISTANCE) {
                            neighborDistances.add(distance);
                        }
                    }
                }                            //
            }
        }
    }

    public ArrayList<PeptideBondDistanceList> neighbors() {
        if (neighbors == null)
            throw new RuntimeException("Neighbors is null.\n The method \" public ArrayList<PeptideBondDistanceList> neighbors(ArrayList<Segment> segments) \" has not been called yet.");
        return neighbors;
    }

    public int numberOfNeighbors() {
        int numberOfNeighbors = 0;
        neighbors();
        for (PeptideBondDistanceList dl : neighbors)
            if (!dl.isEmpty()) numberOfNeighbors++;
        return numberOfNeighbors;
    }
}
