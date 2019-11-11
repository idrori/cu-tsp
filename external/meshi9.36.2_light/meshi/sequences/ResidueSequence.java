/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import java.util.*;

import meshi.util.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

public class ResidueSequence extends MeshiSequence {
    public static final ResidueSequenceCharFilter residueCharFilter = new ResidueSequenceCharFilter();

    public ResidueSequence() {
        super(residueCharFilter);
    }

    public ResidueSequence(String comment) {
        super(comment, residueCharFilter);
    }

    public ResidueSequence(String sequence, String comment) {
        super(residueCharFilter);
        Character weirdChar = weirdChar(sequence);
        if (weirdChar == null) { // that is all characters in the sequence are either valid residue
            // codes, gap or wildcard
            comments.add(comment);
            int number = 1;
            for (int i = 0; i < sequence.length(); i++) {
                SequenceAlignmentCell cell;
                char chr = sequence.charAt(i);
                if (chr == SequenceAlignmentCell.GAP_CHAR)
                    cell = new SequenceAlignmentCell();
                else {
                    cell = new SequenceAlignmentCell(chr, number);
                    number++;
                }
                SequenceAlignmentColumn column = new SequenceAlignmentColumn(cell);
                add(column);
            }
        } else {
            throw new ResidueTypeException(weirdChar, sequence);
        }
    }

    public ResidueSequence(MeshiSequence sequence) {
        this(sequence.comment());
        for (Iterator columns = sequence.iterator(); columns.hasNext();)
            add((SequenceAlignmentColumn) columns.next());
    }

    public ResidueSequence(FastaList fastaList) {
        super(fastaList.get(0), residueCharFilter);
        if (fastaList.size() != 1) throw new RuntimeException("fastaList size = " + fastaList.size());
    }

    public ResidueSequence(FastaList fastaList, int index) {
        super(fastaList.get(index), residueCharFilter);
    }

    public ResidueSequence(ResidueList residueList, String comment) {
        super(comment, residueCharFilter);
        for (Iterator residues = residueList.iterator(); residues.hasNext();) {
            Residue residue = (Residue) residues.next();
            if (!residue.dummy()) {
                char residueName = residue.type().nameOneLetter().charAt(0);
                int residueNumber = residue.number();
                SequenceAlignmentCell cell = new SequenceAlignmentCell(residueName, residueNumber);
                cell.addAttribute(residue);
                SequenceAlignmentColumn column = new SequenceAlignmentColumn(cell);
                add(column);
            }
        }
    }


    private static class ResidueSequenceCharFilter extends SequenceCharFilter {
        public boolean accept(Object obj) {
            Character c = ((Character) obj).charValue();
            if ("XACDEFGHIKLMNOPQRSTVWY-".indexOf(c) >= 0) return true;
            return false;
        }
    }
}
