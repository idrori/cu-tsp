/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.parameters.*;

public class CompositePropensityParameters implements Parameters,
        CompositeTorsionsDefinitions {

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

    public CompositePropensityParameters(int residueType,
                                         SplinedPolynomialsLoader spl) {
        this.residueType = residueType;

        propPolynomial = spl.findPolynomial(residueType, POLYNOMIAL_PHI_PSI_TORSIONS, ALL);
        omniPolynomial = spl.findPolynomial(OMNI, POLYNOMIAL_PHI_PSI_TORSIONS, ALL);

        if (propPolynomial == null) {
            System.out.println("CPP error, unable to find propPolynomial");
            System.exit(1);
        }
        if (omniPolynomial == null) {
            System.out.println("CPP error, unable to find omniPolynomial");
            System.exit(1);
        }
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
            case 2:
                derivVar = PSI;
                break;
        }
        ;

        return propPolynomial.value(derivVar, phi, psi) - omniPolynomial.value(derivVar, phi, psi );
	}
}
