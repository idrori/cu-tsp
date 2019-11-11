/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.ParametersList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.*;

public class RamachandranSidechainParametersList extends ParametersList implements CompositeTorsionsDefinitions {

    /**
     * Loader contains all of the polynomials in the data file.
     */
    private static SplinedPolynomialsLoader spl = null;

    /**
     * Create a new parameters file from the polynomials parameters data file.
     */
    public RamachandranSidechainParametersList(
            String splinedPolynomialsFileName) {
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

        /* create parameters for each amino acid type */
        for (int aac = 0; aac < 20; aac++)
            add(createParameters(aac));
    }

    /**
     * Returns parameters for ResidueTorsions object according to
     * residue type.
     */
    public Parameters parameters(Object Obj) {
        ResidueTorsions residueTorsions = (ResidueTorsions) Obj;

        /* scan internal array got possible residue types */
        for (int i = 0; i < size(); i++) {
            RamachandranSidechainParameters rsp =
                    (RamachandranSidechainParameters) get(i);
            if (rsp.getResidueType() == residueTorsions.getResidueType().ordinal())
                return rsp;
        }

        return null;
    }

    /**
     * Creates a RamachandranSidechainParameters for residue type.
     */
    public Parameters createParameters(int residueType) {
        /* classify residue according to type and call appropriate
           * constructor class. types are grouped according to number
           * of side chain torsion angles.
           */
        if (NUM_SIDECHAIN_TORSIONS[residueType] == 0)
            return new RamachandranSidechainParametersChi0(residueType, spl);
        else if (NUM_SIDECHAIN_TORSIONS[residueType] == 1)
            return new RamachandranSidechainParametersChi1(residueType, spl);
        else if (NUM_SIDECHAIN_TORSIONS[residueType] == 2)
            return new RamachandranSidechainParametersChi2(residueType, spl);
        else if (NUM_SIDECHAIN_TORSIONS[residueType] == 3)
            return new RamachandranSidechainParametersChi3(residueType, spl);
        else if (NUM_SIDECHAIN_TORSIONS[residueType] == 4)
            return new RamachandranSidechainParametersChi4(residueType, spl);
        else
            throw new RuntimeException("unrecognized amino acid type");
    }

    /**
     * Obsolete in RamachandranSidechainParametersList.
     */
    public Parameters createParameters(String line) {
        throw new RuntimeException("RamachandranSidechainEnergy term uses " +
                "SplinedPolynomialsLoader in order to create parameters list" );
	}
}
