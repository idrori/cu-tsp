/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.util.filters.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

import java.util.*;

/**
 * Residue alignment of multiple protein objects, all representing the same protein.
 * Thus, they all have identicle residues up to gaps.
 * The alignment includes only the residues that exists in all models.
 */
public class ResidueAlignmentOfModels extends ResidueAlignment {
    public ResidueAlignmentOfModels(Protein proteins[]) {
        super();
        Iterator[] iterators = new Iterator[proteins.length];

        for (int i = 0; i < proteins.length; i++) {
            comments.add(proteins[i].name());
            iterators[i] = proteins[i].residues().iterator();
        }

        while (allHasNext(iterators)) {
            ResidueAlignmentColumn column = new ResidueAlignmentColumn(proteins.length);
            for (int i = 0; i < proteins.length; i++) {
                column.add(i, new ResidueAlignmentCell((Residue) iterators[i].next()));
            }
            if (!column.hasGap()) add(column);
        }
    }


    private boolean allHasNext(Iterator[] iterators) {
        for (int i = 0; i < iterators.length; i++)
            if (!iterators[i].hasNext()) return false;
        return true;
    }


    public String toString() {
        return (new SequenceAlignment(this)).toString();
    }


}
