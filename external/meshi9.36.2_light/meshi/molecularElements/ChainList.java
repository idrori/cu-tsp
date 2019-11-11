/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.filters.Filter;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A list of residues.
 */
public class ChainList extends ArrayList<Chain> {
    public final Protein protein;

    // --------------------------------- constructors ------------------------------

    /**
     * Empty list.
     */
    public ChainList(Chain chain) {
        this(chain.protein);
        add(chain);
    }
    public ChainList(Protein protein) {
        super();
        this.protein = protein;
    }

    public ChainList(AtomList atomList, ResidueCreator creator, Protein protein) {
        this(protein);
        Atom atom;
        AtomList newList = null;
        String chainName = "XXXX";
        String prevChainName = "XXXX";
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            atom = (Atom) atoms.next();
            if (atom.pdbLine().chain().equals(chainName)) {
                newList.add(atom);
            } else {
                prevChainName = chainName;
                chainName = atom.pdbLine().chain();
                if (newList == null) {
                    newList = new AtomList(atomList.molecularSystem());
                    newList.add(atom);
                } else {
                    add(new Chain(newList, creator, prevChainName, protein));
                    newList = new AtomList(atomList.molecularSystem());
                    newList.add(atom);
                }
            }
        }
        if (newList == null)
            throw new RuntimeException("error 1: Cannot create chains from an empty list of atoms. Protein = " + protein.name());
        if (newList.size() < 1) throw new RuntimeException("error 2: Cannot create chains from an empty list of getAtoms");
        add(new Chain(newList, creator, chainName, protein));
    }


// for Image chains

    public ChainList filter(Filter chainFilter) {
        ChainList newList = new ChainList(protein);
        for (Chain chain : this) {
            if (chainFilter.accept(chain))
                newList.add(chain);
        }
        return newList;
    }


    static class IsChain implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof Chain);
        }
    }


    public String toString() {
        return "chainList with " + size() + " chains";
    }
}
	
