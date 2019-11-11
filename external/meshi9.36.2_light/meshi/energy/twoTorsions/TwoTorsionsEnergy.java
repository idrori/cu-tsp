/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.twoTorsions;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionPair;
import meshi.geometry.TorsionPairList;
import meshi.molecularElements.atoms.Atom;
import meshi.util.filters.Filter;

/**
 * TwoTorsionsEnergy potential - This potential gives the energy values and derivatives for a 2 torsions landscape.
 * <p/>
 * This potential gives a 2D cubic spilne interpolation of the 2 torsions plot.
 * The interpolation is periodic meaning that the value and first derivatives at value (+limit) are equal to those
 * at (-limit), for both torsions.
 * <p/>
 * The parameters of this potential are stored in big files, each relevent to a specific 2-torsion landscape.
 * The names of the files are given in the input as a TwoTorsionsParametersList object.
 * This energy class also gets a list of torsions as input. It choses the ones for which it has parameters
 * (based isOn the torsion names and secondary structures), and construct a list of energy elements.
 */

public class TwoTorsionsEnergy extends SimpleEnergyTerm {
    public static final Filter isTorsionPair = new IsTorsionPair();
    protected TorsionPairList torsionPairList;
    protected DistanceMatrix distanceMatrix;

    public TwoTorsionsEnergy() {
    }

    public boolean evaluatesResidues() {
        return false;
    }

    public TwoTorsionsEnergy(TorsionPairList torsionPairList, DistanceMatrix distanceMatrix,
                             TwoTorsionsParametersList parametersList, EnergyInfoElement info, String comment) {
        super(toArray(distanceMatrix, torsionPairList), parametersList, info);
        this.torsionPairList = torsionPairList;
        this.distanceMatrix = distanceMatrix;
        createElementsList(torsionPairList);
        this.comment = comment;
    }


    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        TorsionPair torsionPair = (TorsionPair) baseElement;
        for (Atom atom : torsionPair.torsion1().atoms())
            if (atom.nowhere()) return null;
        for (Atom atom : torsionPair.torsion2().atoms())
            if (atom.nowhere()) return null;
        return new TwoTorsionsEnergyElement(torsionPair,
                (TwoTorsionsParameters) parameters, weight);
    }

     public void scaleWeight(double factor) {
         super.scaleWeight(factor);
         for (EnergyElement element : elementsList)
             ((TwoTorsionsEnergyElement) element).scaleWeight(factor);
     }
    public void handleMissingParameters(Object obj) {
        throw new RuntimeException("No handle missing parameters at " + comment());
    }

    public static class IsTorsionPair implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof TorsionPair);
        }
    }


}    
