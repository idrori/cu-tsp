/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.molecularElements.atoms.AtomList;
import meshi.sequences.SequenceAlignmentCell;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * A list of residues.
 */
public class ResidueList extends ArrayList<Residue> {

    // --------------------------------- constructors ------------------------------

    private File sourceFile = null;
    public File getSourceFile() {
        return sourceFile;
    }
    public void setSourceFile(File sourceFile){
        this.sourceFile = sourceFile;
    }
    public ResidueList() {
        super();
    }

    public ResidueList(ChainList chains) {
        super();  /*
        for (Iterator chainsIter = chains.iterator(); chainsIter.hasNext();) {
            Chain chain = (Chain) chainsIter.next();
            Iterator residues = chain.iterator();
            Residue residue = (Residue) residues.next();
            if (! residue.dummy())
            throw new RuntimeException("The first residue in each Chain(so as in Residues)  must be dummy!");
            add(residue);

            for (; residues.hasNext();) {
                residue = (Residue) residues.next();
                if (! residue.dummy()) add(residue);
            }
        }

        /*/
        for (Iterator chainsIter = chains.iterator(); chainsIter.hasNext();) {
            Chain chain = (Chain) chainsIter.next();
            for (Iterator residues = chain.iterator(); residues.hasNext();) {
                Residue residue = (Residue) residues.next();
                if (!residue.dummy()) add(residue);
            }
        }     //*/

    }
    //------------------------------------------- methods -----------------------------------
    /**
     * Fetch a residue by its position in the list.
     **/

    /**
     * Fetches a residue by its residue identifier.
     */
    public Residue residue(ResidueIdentifier residueID) {
        Residue residue;
        Iterator residueIterator = iterator();
        while (residueIterator.hasNext()) {
            residue = (Residue) residueIterator.next();
            if (residue.ID().equals(residueID))
                return residue;
        }
        return null;
    }


    /**
     * Extracts the residues accepted by the filter.
     */
    public ResidueList filter(Filter residueFilter) {
        ResidueList newList = new ResidueList();
        for (Residue residue : this) {
            if (residueFilter.accept(residue))
                newList.add(residue);
        }
        return newList;
    }


    public String toString() {
        String out = "";
        Iterator residues = iterator();
        boolean first = true;
        ;
        while (residues.hasNext()) {
            Residue residue = (Residue) residues.next();
            if ((!first) || (!residue.dummy())) {
                if (residue.dummy()) out += SequenceAlignmentCell.WILDCARD_CHAR;
                else out += residue.type().nameOneLetter();
            }
            first = false;
        }
        return out;
    }

    public AtomList atoms() {
        if (size() < 1) throw new RuntimeException("Cannot create an atom list from an empty residue list.");
        return new AtomList(this,get(0).getAtoms().molecularSystem);
    }


    public static class NonDummyFilter implements Filter {
        public boolean accept(Object obj) {
            return !((Residue) obj).dummy();
        }
    }

    public void print() {
        for (Residue r : this)
            System.out.println(r);
    }

    public void print(MeshiWriter writer) {
            for (Residue residue : this)
                if (!residue.dummy())
                    residue.getAtoms().somewhere().print(writer);
        }

    public void sort() {
        Residue[] temp = toArray(new Residue[size()]);
        clear();
        for (Residue r : temp)
            add(r);
    }

    public void freeze() {
        for (Residue residue : this) {
            residue.getAtoms().freeze("In ResidueList");
        }
    }



    /**
     * A secure defrost  method. Does not allow the simulation to continue after getAtoms have melted.
     */
    public void defrost() {
        for (Residue residue : this) {
            residue.getAtoms().defrost();

        }
    }

}
