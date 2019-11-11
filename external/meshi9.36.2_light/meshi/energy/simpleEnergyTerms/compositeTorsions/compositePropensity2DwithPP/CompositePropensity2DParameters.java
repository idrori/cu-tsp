/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity2DwithPP;

import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomial;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.*;

public class CompositePropensity2DParameters implements Parameters, CompositeTorsionsDefinitions {

    /**
     * residue type
     */
    private final int residueType;

    /**
     * propensity polynomial for this residue
     */
    private SplinedPolynomial propPolynomial;

    /**
     * the omni polynomial
     */
    private SplinedPolynomial omniPolynomial;

    public CompositePropensity2DParameters(int residueType,
                                           SplinedPolynomialsLoader spl) {
        this.residueType = residueType;

        propPolynomial = spl.findPolynomial(residueType, POLYNOMIAL_PHI_PSI_TORSIONS, ALL);
        omniPolynomial = spl.findPolynomial(OMNI, POLYNOMIAL_PHI_PSI_TORSIONS, ALL);
    }

    /**
     * returns rotamer's residue type
     */
    public int getResidueType() {
        return residueType;
    }


    public double evaluate(int derivVar, ResidueTorsions resTorsions) {
        double phi = resTorsions.getTorsion(PHI).torsion();
        double psi = resTorsions.getTorsion(PSI).torsion();

        switch (derivVar) {
            case 0:
                break;
            case 1:
                derivVar = PHI;
                break;
        }
        ;

        return propPolynomial.value(derivVar, phi, psi) - omniPolynomial.value(derivVar, phi, psi);
	}
}
