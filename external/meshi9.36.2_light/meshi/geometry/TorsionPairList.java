/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomPairList;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.Filter;

import java.util.ArrayList;

/**
 * A list of torsion pairs, used mainly for the various two torsion energies.
 */
public class TorsionPairList extends ArrayList<TorsionPair> implements Updateable {
    private int numberOfUpdates = 0;

    public TorsionPairList() {
        super();
    }

    /**
     * A TorsionPair list based isOn a Torsion list. Currently, only torsion pairs from the same
     * residue are considered. Mixed torsion pairs from different residues are not treated nor
     * created. Also, only torsions with known biological names (PHI , CHI1, etc.) are treated.
     */

    public TorsionPairList(TorsionList torsions) {
        for (int iTorsion1 = 0; iTorsion1 < torsions.size(); iTorsion1++) {
            Torsion torsion1 = torsions.get(iTorsion1);
            for (int iTorsion2 = iTorsion1 + 1; iTorsion2 < torsions.size(); iTorsion2++) {
                Torsion torsion2 = torsions.get(iTorsion2);
                if ((torsion2.getTorsionResNum() == torsion1.getTorsionResNum()) &&
                        (torsion1.getTorsionCode() > -1) &&
                        (torsion2.getTorsionCode() > -1)) {
                    add(new TorsionPair(torsion1, torsion2));
                    add(new TorsionPair(torsion2, torsion1));
                }
            }
        }
    }


    public void update(int numberOfUpdates) throws UpdateableException {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            int size = size();
            for (int i = 0; i < size; i++) {
                get(i).update(numberOfUpdates);
            }
            this.numberOfUpdates++;
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with TorsionPairList.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }


    public boolean sortable() {
        return false;
    }

    /**
     * Create a torsion pair list from a protein. Note, that many of the torsion pairs created are
     * considered not relevent, such as {PSI , CHI4}.
     */
    public static TorsionPairList createTorsionPairList(Protein protein, DistanceMatrix distanceMatrix) {
        AtomPairList bondList = protein.bonds();
        AngleList angleList = new AngleList(bondList, distanceMatrix);
        TorsionList torsionList = new TorsionList(angleList, distanceMatrix);
        return new TorsionPairList(torsionList);
    }

    public static TorsionPairList createQuickAndDirtyTorsionPairList(Protein protein, DistanceMatrix distanceMatrix) {
        AtomPairList bondList = protein.bonds();
        AngleList angleList = new AngleList(bondList, distanceMatrix);
        TorsionList torsionList = new QuickAndDirtyTorsionList(angleList, distanceMatrix);
        return new TorsionPairList(torsionList);
    }

    public static TorsionPairList createQuickAndDirtyTorsionPairList(Protein protein, DistanceMatrix distanceMatrix, Filter torsionsListFilter) {
        AtomPairList bondList = protein.bonds();
        AngleList angleList = new AngleList(bondList, distanceMatrix);
        TorsionList completeTorsionList = new QuickAndDirtyTorsionList(angleList, distanceMatrix);
        /* filter the list as required */
        TorsionList filteredTorsionList = completeTorsionList.filter(torsionsListFilter);

        return new TorsionPairList(filteredTorsionList);
    }

    public TorsionPairList filter(Filter filter) {
        TorsionPairList out = new TorsionPairList();
        for (TorsionPair tp : this) {
            if (filter.accept(tp)) out.add(tp);
        }
        return out;
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
        for (TorsionPair torsionPair : this)
            torsionPair.setNumberOfUpdates(numberOfUpdates);
    }
}

