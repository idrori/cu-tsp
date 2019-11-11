/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.parameters.*;

public class CompositePropensityParametersList extends ParametersList implements CompositeTorsionsDefinitions {

    private static SplinedPolynomialsLoader spl = null;

    public CompositePropensityParametersList(String splinedPolynomialsFileName) {
        super();

        /* do not load spl if it has already been loaded */
        if (spl == null)
            /* attempt to load polynomials from file */
            try {
                spl = new SplinedPolynomialsLoader(splinedPolynomialsFileName);
            }
            catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException("unable to read polynomials parameters file");
            }

        /* create parameters for each amino acid type and OMNI */
        for (int aac = 0; aac < 20; aac++)
            add(createParameters(aac));
        add(createParameters(OMNI));
    }

    public Parameters createParameters(int residueType) {
        return new CompositePropensityParameters(residueType, spl);
    }

    /**
     * Returns parameters for ResidueTorsions object according to
     * residue type.
     */
    public Parameters parameters(Object Obj) {
        ResidueTorsions residueTorsions = (ResidueTorsions) Obj;

        /* scan internal array got possible residue types */
        for (int i = 0; i < size(); i++) {
            CompositePropensityParameters cpp = (CompositePropensityParameters) get(i);
            if (cpp.getResidueType() == residueTorsions.getResidueType().ordinal())
                return cpp;
        }

        return null;
    }

    public Parameters createParameters(String line) {
        throw new RuntimeException("CompositePropensityEnergy term uses " +
                "SplinedPolynomialsLoader in order to create parameters list");
	}

}
