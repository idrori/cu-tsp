/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomial;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomialsLoader;

/**
 * Ramachandran/Sidechain parameters for amino acids with no
 * sidechain torsion angles: ALA, GLY.
 *
 * @author El-ad David Amir
 */
public class RamachandranSidechainParametersChi0 extends
        RamachandranSidechainParameters {

    public RamachandranSidechainParametersChi0(
            int residueType, SplinedPolynomialsLoader spl) {
        super(residueType);

        /* assemble polynomials for amino acid with no sidechain
           * torsion angle
           */
        polynomials = new SplinedPolynomial[1];
        polynomials[POLYNOMIAL_PHI_PSI] =
                spl.findPolynomial(residueType, POLYNOMIAL_PHI_PSI_TORSIONS, ALL);
    }

    protected boolean legalResidueType() {
        return (getResidueType() == ResidueType.ALA.ordinal() || getResidueType() == ResidueType.GLY.ordinal());
    }

    public double evaluate(int derivVar, ResidueTorsions resTorsions) {
        /* evaluate energy of parameter. generally, can be phrased as:
           * Pr(rot) = Pr(phi,psi)
           *         --------->>>
           * En(rot) = En(phi,psi)
           */
        double phi = resTorsions.getTorsion(PHI).torsion();
        double psi = resTorsions.getTorsion(PSI).torsion();
        SplinedPolynomial enpp =
                polynomials[POLYNOMIAL_PHI_PSI];

        if (derivVar == 0)
            return
                    enpp.value(0, phi, psi);
        else if (derivVar == PHI)
            return
                    enpp.value(1, phi, psi);
        else if (derivVar == PSI)
            return
                    enpp.value(2, phi, psi);
        else
            throw new RuntimeException("derived torsion angle identifier doesn't exist in these parameters");
    }

}
