/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.util.*;
import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;

/**
 * Ramachandran/Sidechain energy element for amino acids with one
 * sidechain torsion angle: CYS, SER, THR, VAL.
 *
 * @author El-ad David Amir
 */
public class RamachandranSidechainEnergyElementChi1
        extends RamachandranSidechainEnergyElement
        implements CompositeTorsionsDefinitions {

    public RamachandranSidechainEnergyElementChi1(
            ResidueTorsions residueTorsions,
            RamachandranSidechainParameters rsp,
            double weight) {
        super(residueTorsions, rsp, weight);
    }

    protected boolean legalResidueType() {
        return (residueTorsions.getResidueType() == ResidueType.CYS ||
                residueTorsions.getResidueType() == ResidueType.SER ||
                residueTorsions.getResidueType() == ResidueType.THR ||
                residueTorsions.getResidueType() == ResidueType.VAL);
    }

    public double evaluate() {
        /* verify energy element is not frozen */
        if (frozen()) return 0.0;
        ResidueTorsionsAttribute rta = (ResidueTorsionsAttribute) residueTorsions.getAttribute(MeshiAttribute.RESIDUE_TORSIONS_ATTRIBUTE);
        if (rta == null) {
            rta = new ResidueTorsionsAttribute();
            residueTorsions.addAttribute(rta);
        }
        /* calcualte energy and derivative */
        double energy = rsp.evaluate(0, residueTorsions);
        rta.phi_deriv = rsp.evaluate(PHI, residueTorsions);
        rta.psi_deriv = rsp.evaluate(PSI, residueTorsions);
        rta.chi_1_deriv = rsp.evaluate(CHI_1, residueTorsions);

        RamachandranSidechainEnergy.sumPerResidueType[residueTypeIndex] += energy;
        RamachandranSidechainEnergy.sm2PerResidueType[residueTypeIndex] += energy * energy;
        residueTorsions.setEnergy(energy);
        /* apply weight */
        energy *= weight;
// 		phi_deriv *= weight;
// 		psi_deriv *= weight;
// 		chi_1_deriv *= weight;

        /* apply force to torsions */
        residueTorsions.applyForce(PHI, -rta.phi_deriv * weight);
        residueTorsions.applyForce(PSI, -rta.psi_deriv * weight);
        residueTorsions.applyForce(CHI_1, -rta.chi_1_deriv * weight);

        monitor(energy, rta.chi_1_deriv);
        return energy;
    }
}
