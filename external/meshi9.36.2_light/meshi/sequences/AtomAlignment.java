/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.parameters.ResidueType;
import meshi.util.Rms;
import meshi.util.filters.Filter;
import meshi.util.filters.KolDichfin;

import java.util.ArrayList;
import java.util.Iterator;

public class AtomAlignment extends ArrayList<AtomAlignmentColumn> {
    public AtomAlignment() {
        super();
    }

    public AtomAlignment(AtomList atomList0, AtomList atomList1) {
        this();
        if (atomList0.size() != atomList1.size()) throw new RuntimeException("List lengths must be equal");
        Iterator atoms1 = atomList1.iterator();
        for (Iterator atoms0 = atomList0.iterator(); atoms0.hasNext();)
            add(new AtomAlignmentColumn((Atom) atoms0.next(), atomList0.comment(),
                    (Atom) atoms1.next(), atomList1.comment()));
    }

//    public AtomAlignment(ResidueAlignment residueAlignment) {
//        this(residueAlignment, new KolDichfin());
//    }

    public AtomAlignment(ResidueAlignment[] residueAlignments, Filter filter) {
        super();
        for (int iAlignment = 0; iAlignment < residueAlignments.length; iAlignment++) {
            if (residueAlignments[iAlignment].size() < 1)
                throw new RuntimeException("Residue alignment of size less than one\n"+residueAlignments[iAlignment]);
            if (residueAlignments[iAlignment].get(0).size() != 2)
                throw new RuntimeException("This constructor supports only pairwise residueAlignments\n"+residueAlignments[iAlignment]);
            for (Iterator columns = residueAlignments[iAlignment].iterator(); columns.hasNext();) {
                ResidueAlignmentColumn column = (ResidueAlignmentColumn) columns.next();
                if (!column.hasGap()) {
                    for (Iterator atoms0 = column.residue0().getAtoms().iterator(); atoms0.hasNext();) {
                        Atom atom0 = (Atom) atoms0.next();
                        if (filter.accept(atom0) & (!atom0.nowhere())) {
                            for (Iterator atoms1 = column.residue1().getAtoms().iterator(); atoms1.hasNext(); ) {
                                Atom atom1 = (Atom) atoms1.next();
                                if ((filter.accept(atom1)) & (atom0.name().equals(atom1.name())) & (!atom1.nowhere()))
                                    add(new AtomAlignmentColumn(atom0, atom1));
                            }
                        }
                    }
                }
            }
        }
    }



    public Rms rms() {
        if (hasGaps()) throw new RuntimeException("Cannot calculate RMS for AtomAlignment with gaps");
        return new Rms(this);
    }

    public Atom atomAt(int coulumnIndex, int rowIndex) {
        return (get(coulumnIndex)).atomAt(rowIndex);
    }

    public AtomList atomList(int row) {
        AtomList out = new AtomList(atomAt(0,row).molecularSystem);
        for (Iterator columns = iterator(); columns.hasNext();) {
            AtomAlignmentColumn column = (AtomAlignmentColumn) columns.next();
            AtomAlignmentCell cell = (AtomAlignmentCell) column.cell(row);
            out.add(cell.atom());
        }
        return out;
    }

    /**
     * Does the alignemnt include gaps? .
     */
    public boolean hasGaps() {
        for (Iterator iter = iterator(); iter.hasNext();) {
            AlignmentColumn column = (AlignmentColumn) iter.next();
            if (column.hasGap()) return true;
        }
        return false;
    }

    public static AtomAlignment alignEquivalentAtoms(AtomList atomListI, AtomList atomListJ) {
        Atom atomI, atomJ;
        String atomNameI, atomNameJ;
        Residue residueI, residueJ;
        ResidueType residueTypeI, residueTypeJ;
        int residueInumber, residueJnumber;

        AtomAlignment out = new AtomAlignment();
        for (int i = 0; i < atomListI.size(); i++) {
            atomI = atomListI.get(i);
            if (!atomI.nowhere()) {
                atomNameI = atomI.name();
                residueI = atomI.residue();
                residueTypeI = residueI.type();
                residueInumber = residueI.number();
                for (int j = 0; j < atomListJ.size(); j++) {
                    atomJ = atomListJ.get(j);
                    if (!atomJ.nowhere()) {
                        atomNameJ = atomJ.name();
                        residueJ = atomJ.residue();
                        residueTypeJ = residueJ.type();
                        residueJnumber = residueJ.number();
                        if ((residueInumber == residueJnumber) &&
                                (residueTypeI == residueTypeJ) &&
                                (atomNameI.equals(atomNameJ))) {
                            out.add(new AtomAlignmentColumn(atomI, atomJ));
                            break;
                        }
                    }
                }
            }
        }
        return out;
    }

    public void print() {
        for (AtomAlignmentColumn column : this)
            System.out.println(column);
    }
}
