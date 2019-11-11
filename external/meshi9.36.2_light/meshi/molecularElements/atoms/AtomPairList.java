/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.atoms;

import meshi.molecularElements.Chain;
import meshi.molecularElements.ChainList;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.parameters.AtomType;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.filters.Filter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

public class AtomPairList extends ArrayList<AtomPair> {
    /**
     * An empty AtomPair list
     */
    public AtomPairList() {
        super();
    }

    public AtomPairList(ChainList chainList, CommandList commands) {
        this();
        Iterator bonds;
        AtomPair bond;
        ResidueList dummies = new ResidueList();
        for (Iterator chains = chainList.iterator(); chains.hasNext();) {
            Chain chain = (Chain) chains.next();
            Residue prevResidue = null;
            boolean nonDummyFound = false;
            for (Residue residue : chain) {
                if (residue.dummy() & nonDummyFound) dummies.add(residue); // this is a chain break
                if (!residue.dummy()) {
                    nonDummyFound = true;
                    add(residue.bonds());
                    if ((prevResidue != null) && (!prevResidue.dummy())) {// we are not in the first residue.
                        if (prevResidue.ID().number() != residue.ID().number() - 1)
                            throw new RuntimeException("Weird consecutive residues " + prevResidue + " and " + residue);

                        if (prevResidue.ID().chain() != residue.ID().chain())
                            throw new RuntimeException("Weird consecutive residues " + prevResidue + " of chain " +
                                    prevResidue.ID().chain() +
                                    " and " + residue + " of chain " + residue.ID().chain());
                        if ((prevResidue.nextAtom() == residue.tail()) &&
                                (residue.prevAtom() == prevResidue.head())) {
                            add(new AtomPair(prevResidue.head(), residue.tail()));
//                            throw new RuntimeException("Just wondered what happens here - Chen "+prevResidue.head()+" "+residue.tail());
                        }
                        else {
                            if (prevResidue.nextAtom() != null)   // Previous residue is already bound.
                                throw new RuntimeException("Binding error 1 - An attempt to bind " + prevResidue + " to " + residue +
                                        " while it is already bound to " + prevResidue.nextAtom());
                            if (residue.prevAtom() != null)   // Current aresidue is already bound.
                                throw new RuntimeException("Binding error 2 - An attempt to bind " + residue + " to " + prevResidue +
                                        " while it is already bound to " + residue.prevAtom());
                            if (prevResidue.head() == null)
                                throw new RuntimeException("Binding error 3 - An attempt to bind " + residue + " to " + prevResidue +
                                        " while the head atom of " + prevResidue + " is null");
                            if (residue.tail() == null)
                                throw new RuntimeException("Binding error 4 - An attempt to bind " + residue +
                                        " with null tail to " + prevResidue);
                            add(prevResidue.head().bond(residue.tail()));
                            prevResidue.setNextAtom(residue.tail());
                            residue.setPrevAtom(prevResidue.head());
                        }
                    }
                }
                prevResidue = residue;
            }
            for (int iAtom = 0; iAtom < chain.atoms().size(); iAtom++) {
                Atom atomI = chain.atoms().get(iAtom);
                if ((atomI.type()== AtomType.CSG) && (!atomI.nowhere())) {
                    for (int jAtom = iAtom+1; jAtom < chain.atoms().size(); jAtom++){
                        Atom atomJ = chain.atoms().get(jAtom);
                        if ((atomJ.type()== AtomType.CSG)&& (!atomJ.nowhere())){
                            double d = atomI.distanceFrom(atomJ);
                            if (d < 2.2) {
                                add(atomI.bond(atomJ));
                                //Utils.println("Disulfide found by geometry\n" + atomI + "\n" + atomJ);
                            }
                        }

                    }
                }
            }
            if (commands != null) {
                if (commands.keyExists("disulfide")) {
                    CommandList disulfideCommands = commands.firstWordFilter("disulfide");
                    for (Command command : disulfideCommands) {
                        int residueNumber1 = command.secondWordInt();
                        int residueNumber2 = command.thirdWordInt();
                        Atom atom1 = null;
                        Atom atom2 = null;
                        for (Atom atom:chain.atoms()) {
                            if ((atom.residueNumber() == residueNumber1) &&
                                (atom.type() == AtomType.CSG)) atom1 = atom;
                            if ((atom.residueNumber() == residueNumber2)  &&
                                (atom.type() == AtomType.CSG)) atom2 = atom;
                        }
                        if ((atom1 == atom2 ) || (atom1 == null))
                            throw new RuntimeException("Weird disulfide command"+command);
                        boolean found = false;
                        for (AtomPair disulfide : this) {
                            if ((disulfide.atom1() == atom1) && (disulfide.atom2() == atom2))
                                found = true;
                            else if ((disulfide.atom1() == atom2) && (disulfide.atom2() == atom1))
                                found = true;
                        }
                        if (!found) {
                            Utils.println("Disulfide found from commands\n"+atom1+"\n"+atom2);
                            add(atom1.bond(atom2));
                        }
                    }
                }
            }
        }
    }


    public boolean add(AtomPair pair) {
        return super.add(pair);
    }
    public boolean add(AtomPairList pairList) {
        boolean out = true;
        for (AtomPair atomPair : pairList)
            out = out & add(atomPair);
        return out;
    }
    public AtomList atomList() {
        if (size() == 0)
            throw new RuntimeException("Cannot create an AtomList from an empty AtomPairsList.");
        AtomList list = new AtomList(get(0).atom1().molecularSystem);
        Atom atom1, atom2;
        for (AtomPair atomPair : this) {
            atom1 = atomPair.atom1();
            atom2 = atomPair.atom2();
            if (!list.contains(atom1)) {
                list.add(atom1);
            }
            if (!list.contains(atom2)) {
                list.add(atom2);
            }
        }
        return list;
    }

    public static class IsAtomPair implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof AtomPair);
        }
    }

    public void printLaconic(String prompt) {
        for (AtomPair atomPair : this) {
            System.out.println(atomPair.laconic(prompt));
        }
    }

    public void sort() {
        AtomPair[] array = toArray(new AtomPair[size()]);
        Arrays.sort(array);
        clear();
        for (AtomPair ap : array)
            add(ap);
    }


    public AtomPairList filter(Filter filter) {
        AtomPairList out = new AtomPairList();
        for (AtomPair ap : this) {
            if (filter.accept(ap)) out.add(ap);
        }
        return out;
    }

}
    