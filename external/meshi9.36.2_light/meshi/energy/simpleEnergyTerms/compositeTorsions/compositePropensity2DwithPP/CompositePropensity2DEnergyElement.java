/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity2DwithPP;

import meshi.energy.EnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.parameters.*;

public class CompositePropensity2DEnergyElement
        extends EnergyElement
        implements CompositeTorsionsDefinitions {

    public final ResidueTorsions residueTorsions;
    public final CompositePropensity2DParameters cpp;
    protected double weight;

    public CompositePropensity2DEnergyElement(
            ResidueTorsions residueTorsions,
            CompositePropensity2DParameters cpp,
            double weight) {
        this.residueTorsions = residueTorsions;
        this.cpp = cpp;
        this.weight = weight;

        setAtoms();
        updateFrozen();
    }

    protected void setAtoms() {
        int[] interestingTorsions = {PHI, PSI};
        atoms = residueTorsions.getAtoms(interestingTorsions);
    }

    public double evaluate() {
        /* verify energy element is not frozen */
        if (frozen()) return 0.0;

        /* calcualte energy and derivative */
        double energy = cpp.evaluate(0, residueTorsions);
        double phi_deriv = cpp.evaluate(PHI, residueTorsions);
        double psi_deriv = cpp.evaluate(PSI, residueTorsions);

        /* apply weight */
        energy *= weight;
        phi_deriv *= weight;
        psi_deriv *= weight;

        /* apply force to torsions */
        residueTorsions.applyForce(PHI, -phi_deriv);
        residueTorsions.applyForce(PSI, -psi_deriv);

        return energy;
    }
}
