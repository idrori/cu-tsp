/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.util.MeshiAttribute;

public class CompositePropensityEnergyElement
        extends EnergyElement
        implements CompositeTorsionsDefinitions {

    protected ResidueTorsions residueTorsions;
    protected CompositePropensityParameters cpp;
    protected double weight;
    int residueTypeIndex;

    public CompositePropensityEnergyElement(
            ResidueTorsions residueTorsions,
            CompositePropensityParameters cpp,
            double weight) {
        this.residueTorsions = residueTorsions;
        this.cpp = cpp;
        this.weight = weight;
        setAtoms();
        updateFrozen();
        residueTypeIndex = residueTorsions.getResidueType().ordinal();
        CompositePropensityEnergy.numberOfResiduesPerType[residueTypeIndex]++;
    }

    protected void setAtoms() {
        atoms = residueTorsions.getAtoms();
    }

    public double evaluate() {
        /* verify energy element is not frozen */
        if (frozen()) return 0.0;
        ResidueTorsionsPropensityAttribute rta = (ResidueTorsionsPropensityAttribute) residueTorsions.getAttribute(MeshiAttribute.RESIDUE_TORSIONS_ATTRIBUTE);
        if (rta == null) {
            rta = new ResidueTorsionsPropensityAttribute();
            residueTorsions.addAttribute(rta);
        }

        /* calcualte energy and derivative */
        double energy = cpp.evaluate(0, residueTorsions);
        rta.phi_deriv = cpp.evaluate(PHI, residueTorsions);
        rta.psi_deriv = cpp.evaluate(PSI, residueTorsions);

        CompositePropensityEnergy.sumPerResidueType[residueTypeIndex] += energy;
        CompositePropensityEnergy.sm2PerResidueType[residueTypeIndex] += energy * energy;

        residueTorsions.setEnergy(energy);
        /* apply weight */
        energy *= weight;

        /* apply force to torsions */
        residueTorsions.applyForce(PHI, -rta.phi_deriv * weight);
        residueTorsions.applyForce(PSI, -rta.psi_deriv * weight);


        return energy;
    }
}
