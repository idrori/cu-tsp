/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.bond;

import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfoElement;

import java.util.*;

/**
 * Bond energy term.
 * Has the general form <b> Eb = SIGMAi(Ki(Di-D0i)^2)</b>
 * where <b>Di</b>
 * is the distance between the two bonded getAtoms, <b> D0i </b>
 * is their expected average distance (depends isOn their types) and <b>Ki</b>
 * is a force constant that again, depends isOn the atom types.<br>
 * This class is used for both calculating the bond-energy term of an energy function
 * and for updating the forces isOn each atom accordingly.<b>
 * It is assumed that the list of bonds is constant during the simulation. That is
 * no bonds are made or broken.
 */
public class BondEnergy extends SimpleEnergyTerm {
    /**
     * The constructor associates any bond with its parameters.
     */
    protected DistanceMatrix distanceMatrix;
    private  boolean simpleFlag;
    public BondEnergy() {
    }

    public BondEnergy(AtomPairList bondList,
                      DistanceMatrix distanceMatrix,
                      BondParametersList parametersList,
                      EnergyInfoElement info, boolean simpleFlag) {
        super(toArray(distanceMatrix), parametersList, info);
        comment = "Bond";
        this.distanceMatrix = distanceMatrix;
        this.simpleFlag = simpleFlag;
        createElementsList(bondList);
    }

    public double weight() {
        return weight;
    }

    public void setWeight(double newWeight) {
        weight = newWeight;
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        AtomPair atomPair = (AtomPair) baseElement;
        if (atomPair.atom1().nowhere() || atomPair.atom2().nowhere()) return null;
        if ((distanceMatrix.distance(atomPair.atom1(), atomPair.atom2()) == null) &&
                (distanceMatrix.distance(atomPair.atom2(), atomPair.atom1()) == null)) return null;
        return new BondEnergyElement(atomPair, parameters, distanceMatrix, weight,simpleFlag);
    }

    public void scaleWeight(double factor) {
        super.scaleWeight(factor);
        for (EnergyElement element : elementsList)
            ((BondEnergyElement) element).scaleWeight(factor);
    }

    public void removeBadBonds(double farAway) {
        for (Iterator elements = elementsList.iterator(); elements.hasNext();) {
            BondEnergyElement element = (BondEnergyElement) elements.next();
            Atom atom1 = element.atom1();
            Atom atom2 = element.atom2();
            if (atom1.distanceFrom(atom2) > farAway)
                element.turnOff();
            else element.turnOn();
        }
    }
}
	
	
