/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.EnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.molecularElements.*;
import meshi.parameters.*;

/**
 * Encapsulation of RamachandranSidechain energy value for a single residue.
 * Much like the parameters for this energy function, each residue type
 * has its own class of energy element.
 *
 * @author El-ad David Amir
 */
public abstract class RamachandranSidechainEnergyElement extends EnergyElement {

    protected ResidueTorsions residueTorsions;
    protected RamachandranSidechainParameters rsp;
    protected double weight;

    public double weight() {
        return weight;
    }

    protected int residueTypeIndex;

    public RamachandranSidechainEnergyElement(
            ResidueTorsions residueTorsions,
            RamachandranSidechainParameters rsp,
            double weight) {
        this.residueTorsions = residueTorsions;
        this.rsp = rsp;
        this.weight = weight;

        if (!legalResidueType())
            throw new RuntimeException("energy element residue type mismatch");

        setAtoms();
        updateFrozen();
        residueTypeIndex = residueTorsions.getResidueType().ordinal();
        RamachandranSidechainEnergy.numberOfResiduesPerType[residueTypeIndex]++;
    }

    public Residue residue() {
        return residueTorsions.getResidue();
    }

    public String toString() {
        return "RamachandranSidechainEnergyElement element of " + residueTorsions.toString();
    }

    public void scaleWeight(double factor) {
        weight *= factor;
    }

    protected void setAtoms() {
        atoms = residueTorsions.getAtoms();
    }

    /**
     * Reports energy values. Currently switched off.
     */
    protected void monitor(double energy, double... derivs) {
        if (false) {
            System.out.println("[[BEGIN]] RamachandranSidechainEnergy monitor information");
            System.out.println(residueTorsions);
            System.out.println("total energy = " + energy);
            for (int i = 0; i < derivs.length; i++)
                System.out.println("deriv #" + (i + 1) + " = " + derivs[i]);
            System.out.println("[[END]] RamachandranSidechainEnergy monitor information");
        }
    }

    /**
     * verifies residue type is a legal residue types for class.
     */
    protected abstract boolean legalResidueType();
	
}
