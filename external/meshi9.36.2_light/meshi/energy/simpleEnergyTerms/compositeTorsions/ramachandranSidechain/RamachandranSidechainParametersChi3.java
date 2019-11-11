/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomial;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomialsLoader;

/**
 * Ramachandran/Sidechain parameters for amino acids with three
 * sidechain torsion angles: GLU, MET, GLN.
 *
 * @author El-ad David Amir
 */
public class RamachandranSidechainParametersChi3 extends
        RamachandranSidechainParameters {

    public RamachandranSidechainParametersChi3(
            int residueType, SplinedPolynomialsLoader spl) {
        super(residueType);

        /* assemble polynomials for amino acid with three sidechain
           * torsion angles
           */
        polynomials = new SplinedPolynomial[5];
        polynomials[POLYNOMIAL_PHI_PSI_CHI_1] =
                spl.findPolynomial(residueType, POLYNOMIAL_PHI_PSI_CHI_1_TORSIONS, ALL);
        polynomials[POLYNOMIAL_CHI_1_CHI_2] =
                spl.findPolynomial(residueType, POLYNOMIAL_CHI_1_CHI_2_TORSIONS, ALL);
        polynomials[POLYNOMIAL_CHI_1] =
                spl.findPolynomial(residueType, POLYNOMIAL_CHI_1_TORSIONS, ALL);
        polynomials[POLYNOMIAL_CHI_1_CHI_3] =
                spl.findPolynomial(residueType, POLYNOMIAL_CHI_1_CHI_3_TORSIONS, ALL);
    }

    protected boolean legalResidueType() {
        return (getResidueType() == ResidueType.GLU.ordinal() || getResidueType() == ResidueType.MET.ordinal() ||
                getResidueType() == ResidueType.GLN.ordinal());
    }

    public double evaluate(int derivVar, ResidueTorsions resTorsions) {
        /* evaluate energy of parameter. generally, can be phrased as:
           * Pr(rot) = Pr(phi,psi,chi_1) * Pr(chi_2|chi_1) * Pr(chi_3|chi_1)
           *         = (Pr(phi,psi,chi_1) *
           *           (Pr(chi_1,chi_2) / Pr(chi_1)) *
           *           (Pr(chi_1,chi_3) / Pr(chi_1))
           *         --------->>>
           * En(rot) = En(phi,psi,chi_1) +
           *           En(chi_1,chi_2) - En(chi_1) +
           *           En(chi_1,chi_3) - En(chi_1)
           */
        double phi = resTorsions.getTorsion(PHI).torsion();
        double psi = resTorsions.getTorsion(PSI).torsion();
        double chi_1 = resTorsions.getTorsion(CHI_1).torsion();
        double chi_2 = resTorsions.getTorsion(CHI_2).torsion();
        double chi_3 = resTorsions.getTorsion(CHI_3).torsion();
        SplinedPolynomial enpp1 =
                polynomials[POLYNOMIAL_PHI_PSI_CHI_1];
        SplinedPolynomial en12 =
                polynomials[POLYNOMIAL_CHI_1_CHI_2];
        SplinedPolynomial en1 =
                polynomials[POLYNOMIAL_CHI_1];
        SplinedPolynomial en13 =
                polynomials[POLYNOMIAL_CHI_1_CHI_3];

        if (derivVar == 0)
            return
                    enpp1.value(0, phi, psi, chi_1) +
                            en12.value(0, chi_1, chi_2) -
                            en1.value(0, chi_1) +
                            en13.value(0, chi_1, chi_3) -
                            en1.value(0, chi_1);
        else if (derivVar == PHI)
            return
                    enpp1.value(1, phi, psi, chi_1);
        else if (derivVar == PSI)
            return
                    enpp1.value(2, phi, psi, chi_1);
        else if (derivVar == CHI_1)
            return
                    enpp1.value(3, phi, psi, chi_1) +
                            en12.value(1, chi_1, chi_2) -
                            en1.value(1, chi_1) +
                            en13.value(1, chi_1, chi_3) -
                            en1.value(1, chi_1);
        else if (derivVar == CHI_2)
            return
                    en12.value(2, chi_1, chi_2);
        else if (derivVar == CHI_3)
            return
                    en13.value(2, chi_1, chi_3);
        else
            throw new RuntimeException("derived torsion angle identifier doesn't exist in these parameters");
    }

}
