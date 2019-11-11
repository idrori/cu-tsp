/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.parameters.*;
import meshi.util.Utils;

public class RamachandranParametersList extends ParametersList
        implements CompositeTorsionsDefinitions {

    private SplinedPolynomialsLoader spl = null;

    public RamachandranParametersList(String splinedPolynomialsFileName) {
        super();

        /* attempt to load polynomials from file */
        try {
            Utils.println("Loading " + this + " parameters from " + splinedPolynomialsFileName);
            spl = new SplinedPolynomialsLoader(splinedPolynomialsFileName);
        }
        catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("unable to read polynomials parameters file");
        }

        /* create parameters for each amino acid type and OMNI */
        for (int aac = 0; aac <= OMNI; aac++)
            add(createParameters(aac));
    }

    public Parameters createParameters(int residueType) {
        return new RamachandranParameters(residueType, spl);
    }

    /**
     * Returns parameters for ResidueTorsions object according to
     * residue type.
     */
    public Parameters parameters(Object Obj) {
        ResidueTorsions residueTorsions = (ResidueTorsions) Obj;

        /* scan internal array got possible residue types */
        for (int i = 0; i < size(); i++) {
            RamachandranParameters cpp = (RamachandranParameters) get(i);
            if (cpp.getResidueType() == residueTorsions.getResidueType().ordinal())
                return cpp;
        }

        return null;
    }

    public Parameters createParameters(String line) {
        throw new RuntimeException("RamachandranEnergy term uses " +
                "SplinedPolynomialsLoader in order to create parameters list");
	}

}
