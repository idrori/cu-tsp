/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomial;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomialsLoader;

/**
 * Ramachandran/Sidechain parameters for amino acids with one
 * sidechain torsion angle: CYS, SER, THR, VAL.
 *
 * @author El-ad David Amir
 */
public class RamachandranSidechainParametersChi1 extends
        RamachandranSidechainParameters {

    public RamachandranSidechainParametersChi1(
            int residueType, SplinedPolynomialsLoader spl) {
        super(residueType);

        /* assemble polynomials for amino acid with one sidechain
           * torsion angle
           */
        polynomials = new SplinedPolynomial[2];
        polynomials[POLYNOMIAL_PHI_PSI_CHI_1] =
                spl.findPolynomial(residueType, POLYNOMIAL_PHI_PSI_CHI_1_TORSIONS, ALL);
    }

    protected boolean legalResidueType() {
        return (getResidueType() == ResidueType.CYS.ordinal() || getResidueType() == ResidueType.SER.ordinal() ||
                getResidueType() == ResidueType.THR.ordinal() || getResidueType() == ResidueType.VAL.ordinal());
    }

    public double evaluate(int derivVar, ResidueTorsions resTorsions) {
        /* evaluate energy of parameter. generally, can be phrased as:
           * Pr(rot) = Pr(phi,psi,chi_1)
           *         --------->>>
           * En(rot) = En(phi,psi,chi_1)
           */
        double phi = resTorsions.getTorsion(PHI).torsion();
        double psi = resTorsions.getTorsion(PSI).torsion();
        double chi_1 = resTorsions.getTorsion(CHI_1).torsion();
        SplinedPolynomial enpp1 =
                polynomials[POLYNOMIAL_PHI_PSI_CHI_1];

        if (derivVar == 0)
            return
                    enpp1.value(0, phi, psi, chi_1);
        else if (derivVar == PHI)
            return
                    enpp1.value(1, phi, psi, chi_1);
        else if (derivVar == PSI)
            return
                    enpp1.value(2, phi, psi, chi_1);
        else if (derivVar == CHI_1)
            return
                    enpp1.value(3, phi, psi, chi_1);
        else
            throw new RuntimeException("derived torsion angle identifier doesn't exist in these parameters");
    }

}
