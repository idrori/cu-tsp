/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran;

import meshi.energy.EnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;


import meshi.util.Utils;

public class RamachandranEnergyElement
        extends EnergyElement
        implements CompositeTorsionsDefinitions {

    private double sourceCosPhi = -99999, sourceCosPsi = -99999;
    public static final double THRESHOLD = 6.18;
    public static final double FACTOR = 1000.0;
    private boolean strict;
    public final ResidueTorsions residueTorsions;
    public final RamachandranParameters cpp;
    protected double weight;

    public double weight() {
        return weight;
    }

    public RamachandranEnergyElement(
            ResidueTorsions residueTorsions,
            RamachandranParameters cpp,
            double weight) {
        this.residueTorsions = residueTorsions;
        this.cpp = cpp;
        this.weight = weight;

        setAtoms();
        updateFrozen();
        strict = false;
    }

    public String toString() {
        return "RamachandranEnergyElement element of " + residueTorsions.toString();
    }

    protected void setAtoms() {
        int[] interestingTorsions = {PHI, PSI};
        atoms = residueTorsions.getAtoms(interestingTorsions);
    }

    public double evaluate() {
        double cosPhi, dCosPhi, cosPsi, dCosPsi, expPhi, expPsi, sinPhi, sinPsi;
        /* verify energy element is not frozen */
        if (frozen()) return 0.0;
        double energy, phi_deriv, psi_deriv;
        /* calculate energy and derivative */
        energy = cpp.evaluate(0, residueTorsions);
        phi_deriv = cpp.evaluate(PHI, residueTorsions);
        psi_deriv = cpp.evaluate(PSI, residueTorsions);
        if (strict) {
            cosPhi = Math.cos(residueTorsions.phi().torsion());
            cosPsi = Math.cos(residueTorsions.psi().torsion());
            sinPhi = Math.sin(residueTorsions.phi().torsion());
            sinPsi = Math.sin(residueTorsions.psi().torsion());
            dCosPhi = sourceCosPhi - cosPhi;
            dCosPsi = sourceCosPsi - cosPsi;
            expPhi = Math.exp(-dCosPhi * dCosPhi * 10);
            expPsi = Math.exp(-dCosPsi * dCosPsi * 10);
            energy += FACTOR * (expPhi + expPsi);
            phi_deriv += FACTOR * expPhi * -20 * dCosPhi * sinPhi;
            psi_deriv += FACTOR * expPsi * -20 * dCosPsi * sinPsi;
        }
        energy *= weight;
        phi_deriv *= weight;
        psi_deriv *= weight;

        /* apply force to torsions */
        residueTorsions.applyForce(PHI, -phi_deriv);
        residueTorsions.applyForce(PSI, -psi_deriv);

        return energy;
    }

    public void reset() {
        strict = false;
        double e = evaluate() / weight;
        if (e > THRESHOLD) {
            strict = true;
            Utils.println(this + " is strict");
        } else strict = false;
        sourceCosPhi = Math.cos(residueTorsions.phi().torsion());
        sourceCosPsi = Math.cos(residueTorsions.psi().torsion());
    }
}
