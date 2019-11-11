/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran;

import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.parameters.*;

public class RamachandranParameters implements Parameters, CompositeTorsionsDefinitions {

    /**
     * residue type
     */
    private final int residueType;

    /**
     * polynomial for this residue
     */
    private SplinedPolynomial propPolynomial;

    public RamachandranParameters(int residueType,
                                  SplinedPolynomialsLoader spl) {
        this.residueType = residueType;

        propPolynomial = spl.findPolynomial(residueType, POLYNOMIAL_PHI_PSI_TORSIONS, ALL);
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

        return propPolynomial.value(derivVar, phi, psi);
	}
}
