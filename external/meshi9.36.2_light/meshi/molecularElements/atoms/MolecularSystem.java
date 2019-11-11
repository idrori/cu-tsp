/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.atoms;

import meshi.energy.TotalEnergy;
import meshi.geometry.Coordinates;
import meshi.geometry.DistanceMatrix;
import meshi.parameters.AtomType;
import meshi.util.Terminator;
import meshi.util.Utils;

import java.util.ArrayList;

public class MolecularSystem extends ArrayList<AtomCore> {
    public final int ID;
    private static int numberOfMolecularSystems = 0;
    private int numberOfAtoms = 0;
    private TotalEnergy currentTotalEnergy;
    private DistanceMatrix distanceMatrix;
    private final DistanceMatrix.DistanceMatrixType distanceMatrixType;
    private Terminator terminator;
    private static boolean debug = true;
    //private static TotalEnergy currentTotalEnergy;

    public MolecularSystem(DistanceMatrix.DistanceMatrixType distanceMatrixType) {
        super();
        ID = numberOfMolecularSystems;
        numberOfMolecularSystems++;
        currentTotalEnergy = null;
        distanceMatrix = null;
        this.distanceMatrixType = distanceMatrixType;
        terminator = new Terminator();
    }
    public MolecularSystem() {
        this(DistanceMatrix.DistanceMatrixType.STANDARD);
    }

    public Terminator terminator() {return terminator;}
    public void register(TotalEnergy totalEnergy) {
        if (this.currentTotalEnergy == null) this.currentTotalEnergy = totalEnergy;
        else {
            this.currentTotalEnergy.off();
            this.currentTotalEnergy = totalEnergy;
        }
        getDistanceMatrix();
    }


    public DistanceMatrix getDistanceMatrix() {
        return getDistanceMatrix(DistanceMatrix.DistanceMatrixType.STANDARD);
    }

    public DistanceMatrix getDistanceMatrix(DistanceMatrix.DistanceMatrixType distanceMatrixType) {
        if (distanceMatrix == null) createDistanceMatrix("The first distance matrix for this molecular system.",distanceMatrixType);
        else
            if (distanceMatrix.type != distanceMatrixType) {
                createDistanceMatrix("A new DistanceMatrix of type"+distanceMatrixType, distanceMatrixType);
            }

        return distanceMatrix;
    }

    public void createDistanceMatrix() {
        createDistanceMatrix(null);
    }

    public void createDistanceMatrix(String message) {
        createDistanceMatrix(message,DistanceMatrix.DistanceMatrixType.STANDARD);
    }

    public void createDistanceMatrix(String message,DistanceMatrix.DistanceMatrixType distanceMatrixType) {
        if (debug && (message != null))
            Utils.println(this+" Creating new DistanceMatrix because: " + message);

        terminator.kill(message+"\nA new DistanceMatrix was created.");
        terminator = new Terminator();
                  switch (distanceMatrixType) {
                    case STANDARD: distanceMatrix = new DistanceMatrix(this);
                                   break;
                    case LONG_DISTANCES: distanceMatrix = new DistanceMatrix(this,distanceMatrixType);
                                         break;
                    default:throw new RuntimeException("do something "+distanceMatrixType);
                }
    }

    public TotalEnergy currentTotalEnergy() {return currentTotalEnergy;}

    /*public static void setCurrentMolecularSystem(MolecularSystem ms) {
        currentMolecularSystem = ms;
    } */


    public String toString() {
        return "Molecular System " + ID + "     ";
    }

    protected AtomCore createAtomCore(Atom atom, AtomType type, AtomStatus status, double x, double y, double z) {
        if (currentTotalEnergy != null) throw new RuntimeException("Cannot create more getAtoms in this MolecularSystem an energy function already uses it.");
        if (distanceMatrix != null) throw new RuntimeException("Cannot create more getAtoms in this MolecularSystem a distance matrix already uses it.");
        AtomCore newAtomCore = new AtomCore(atom, type, status, size(), x, y, z);
        if ((x<= Coordinates.NOWHERE_CONST) && (status != AtomStatus.NOWHERE))
            throw new RuntimeException("This is weird\n"+x+" "+newAtomCore);
        add(newAtomCore);
        terminator.kill("An atom was added.");
        return newAtomCore;
    }

    public AtomList toAtomList() {
        AtomList out = new AtomList(this);
        for (AtomCore atomCore : this)
            out.add(atomCore.atom);
        return out;
    }


}
