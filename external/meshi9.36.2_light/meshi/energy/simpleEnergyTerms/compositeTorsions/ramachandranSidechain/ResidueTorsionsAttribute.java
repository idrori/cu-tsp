/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.util.*;

public class ResidueTorsionsAttribute implements MeshiAttribute {
    double phi_deriv = 0;
    double psi_deriv = 0;
    double chi_1_deriv = 0;
    double chi_2_deriv = 0;
    double chi_3_deriv = 0;
    double chi_4_deriv = 0;

    public int key() {
        return RESIDUE_TORSIONS_ATTRIBUTE;
    }
}
