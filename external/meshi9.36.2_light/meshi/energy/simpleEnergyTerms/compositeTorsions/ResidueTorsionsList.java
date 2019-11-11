/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions;

import meshi.energy.simpleEnergyTerms.*;


import meshi.geometry.*;
import meshi.molecularElements.*;
import meshi.util.*;
import meshi.parameters.*;
import meshi.util.filters.*;

import java.util.*;

/**
 * An updateable list of ResidueTorsions.
 *
 * @author El-ad David Amir
 */
public class ResidueTorsionsList extends ArrayList<ResidueTorsions> implements Updateable {

    /* updates counter for compatability with Updateable */
    int numberOfUpdates;


    public ResidueTorsionsList() {
        super();
    }

    /**
     * Creates ResidueTorsionsList from protein.
     */
    public ResidueTorsionsList(Protein protein, DistanceMatrix dm) {
        this();

        /* create list of all torsions in protein */
        TorsionList torsionList =
                TorsionList.createQuickAndDirtyTorsionList(protein, dm);

        /* create all residue torsions and add them to list */
        for (Residue res : protein.residues()) {
            if (res.dummy())
                continue;

            ResidueTorsions resTorsions = new ResidueTorsions(res, torsionList);

            /* verify the protein has all the torsion angles required
                * by this residue.
                */
            if (resTorsions.hasAllTorsions())
                add(resTorsions);
        }
    }

    public void update(int numberOfUpdates) throws UpdateableException {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            for (int i = 0; i < size(); i++)
                ((ResidueTorsions) get(i)).update(numberOfUpdates);
            this.numberOfUpdates++;
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("incorrect update number");
    }

    /**
     * A method to filter residues that have {PHI,PSI} torsions.
     * This method removes residues at chain termini or Ca traces, but keeps
     * incomplete residues.
     */
    public ResidueTorsionsList filterPhiPsiResidues() {
        Filter filter = new PhiPsiResidueFilter();
        ResidueTorsionsList out = new ResidueTorsionsList();
        for (ResidueTorsions rt : this)
            if (filter.accept(rt)) out.add(rt);
        return out;
    }

    public static class PhiPsiResidueFilter implements Filter, CompositeTorsionsDefinitions {
        public boolean accept(Object obj) {
            ResidueTorsions tors = (ResidueTorsions) obj;
            return ((tors.getTorsion(PHI) != null) &&
                    (tors.getTorsion(PSI) != null));
        }
    }

    /**
     * A method to filter for incomplete residues. Note that the resulting list
     * CONTAINS all the INCOMPLETE residues of the protein. It is inverse of the
     * filterCompleteResidues() method.
     */
    public ResidueTorsionsList filterIncompleteResidues() {
        return (ResidueTorsionsList) Utils.filter(this, new IncompleteResidueFilter(), new ResidueTorsionsList());
    }

    public static class IncompleteResidueFilter implements Filter {
        public boolean accept(Object obj) {
            ResidueTorsions tors = (ResidueTorsions) obj;
            return !tors.hasAllTorsions();
        }
    }

    /**
     * A method to filter residues that have all torsions
     */
    public ResidueTorsionsList filterCompleteResidues() {
        return (ResidueTorsionsList) Utils.filter(this, new CompleteResidueFilter(), new ResidueTorsionsList());
    }

    public static class CompleteResidueFilter implements Filter {
        public boolean accept(Object obj) {
            ResidueTorsions tors = (ResidueTorsions) obj;
            return tors.hasAllTorsions();
        }
    }

    /**
     * A method to tag as Pre-Proline all the residues that are such.
     */
    private void tagPreProline() {
        for (int c = 0; c < size(); c++) {
            ResidueTorsions thisRT = (ResidueTorsions) get(c);
            int resnum = thisRT.getResidueNumber();
            ResidueTorsions nextRT = findResidueInList(resnum + 1);
            if ((nextRT != null) && (nextRT.getResidueType() == ResidueType.PRO) && //It's a pre-Pro
                    (thisRT.getResidueType() != ResidueType.PRO) && (thisRT.getResidueType() != ResidueType.GLY)) //This is not gly or pro
                thisRT.setToPP();
        }
    }

    public ResidueTorsions findResidueInList(int num) {
        for (int c = 0; c < size(); c++)
            if (((ResidueTorsions) get(c)).getResidueNumber() == num)
                return (ResidueTorsions) get(c);
        return null;
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
        for (ResidueTorsions residueTorsions : this) {
            residueTorsions.setNumberOfUpdates(numberOfUpdates);
}
    }
}
