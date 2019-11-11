/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.molecularElements.Chain;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.aligner.AlignmentScheme;
import meshi.sequences.aligner.DpMatrix;
import meshi.sequences.aligner.SubstitutionScorer;
import meshi.util.MeshiAttribute;

import java.util.Iterator;

/**
 * Sequence is a single row alignment.
 */
public class MeshiSequence extends SequenceAlignment {
    public final SequenceCharFilter charFilter;
    public static final char UNKNOWN = '_';

    public MeshiSequence(String comment) {
        this(new kolDichfin());
        comments.add(comment);
    }

    public MeshiSequence(SequenceCharFilter charFilter) {
        super();
        this.charFilter = charFilter;
    }

    public MeshiSequence(MeshiSequence source) {
        super();
        charFilter = source.charFilter;
        comments.add(source.comments.get(0));
        for (SequenceAlignmentColumn sourceColumn : source) {
            add(new SequenceAlignmentColumn(sourceColumn));
        }
    }

    public MeshiSequence(String comment, SequenceCharFilter charFilter) {
        this(charFilter);
        comments.add(comment);
    }

    public MeshiSequence(Chain chain, String comment) {
        this(getChainSequence(chain), comment, new kolDichfin());
    }

    public MeshiSequence(String sequence, String comment) {//INTERESTING Tommer 24.9.14
        this(sequence, comment, new kolDichfin());
    }


    public MeshiSequence(String sequence, String comment, SequenceCharFilter charFilter) {
        this(charFilter);
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
        }
    }
    public MeshiSequence(String[] temp, String comment, SequenceCharFilter charFilter) {
        this(charFilter);
        String sequence = "";
        for (String s:temp)
            sequence +=  s.charAt(0);
        Character weirdChar = weirdChar(sequence);
        if (weirdChar != null)  // that is not all characters in the sequence are either valid residue codes, gap or wildcard
            throw new RuntimeException("This is weird.") ;
        comments.add(comment);
        sequence = "";
        for (String s:temp)   {
            sequence +=  s.charAt(0);
            int number = 1;
            SequenceAlignmentCell cell;
            char chr = s.charAt(0);
            if (chr == SequenceAlignmentCell.GAP_CHAR)
                    cell = new SequenceAlignmentCell();
            else {
                    cell = new SequenceAlignmentCell(chr, number,s.substring(1));
                    number++;
                }
            SequenceAlignmentColumn column = new SequenceAlignmentColumn(cell);
            add(column);
        }
    }


    public boolean add(SequenceAlignmentColumn column) {
        char c = column.getChar(0);
        Character cc = new Character(c);
        if (!charFilter.accept(cc))
            throw new RuntimeException("Cannot add " + column.getChar(0) + " to " + getClass());
        return super.add(column);
    }

    public Character weirdChar(String sequence) {
        for (int i = 0; i < sequence.length(); i++) {
            char c = sequence.charAt(i);
            if ((c != SequenceAlignmentCell.GAP_CHAR) &
                    (c != SequenceAlignmentCell.WILDCARD_CHAR)) {
                Character cc = new Character(c);
                if (!charFilter.accept(cc)) return cc;
            }
        }
        return null;
    }

    public String comment() {
        if (comments.size() > 0) return comments.get(0);
        throw new RuntimeException("No comment to this Sequence object - " +
                "probably an indication of a problem.");
    }


    public MeshiSequence(MeshiSequence source, SequenceCharFilter charFilter) {
        this(source.comment(), charFilter);
        for (SequenceAlignmentColumn sac : source)
            add(sac);
    }

    public int startsIn() {
        int i;
        boolean found = false;
        for (i = 0; i < size() && (!found); i++){//modified Tommer 11.9.14
        	//System.out.print("\n (startsIn): "+i);
            if (!cell(i).gap()) found = true;
        }
        return (cell(i-1).number - i);//modified by Tommer 11.9.14!! BUG KICKED?
    }

    public String toString() {//Amended for checking by Tommer 1.9.14
        if (comments.size() == 0) {
            if (size() == 0)
                throw new RuntimeException("Weird empty sequence without comment");
            throw new RuntimeException("Weird sequence without comment");
        }
        int commentEnd = comment().indexOf("starts in ");
        if (commentEnd == -1) commentEnd = comment().length();
        String comment = comment().substring(0, commentEnd) + "; starts in " + startsIn();
        String out="";//changed 1.9.14
        //if (comment.startsWith(">")) out = comment;
        //else out = "> " + comment;
        for (int i = 0; i < size(); i++) {
            //if (i % 50 == 0) out += "\n";
            out += charAt(i);
        }
        out += "*" + "\n";
        
        return out;
    }


    public SequenceAlignmentCell cell(int index) {
    	//System.out.print("  (cell):"+ index+" ");
        if (index > size() - 1)
            throw new RuntimeException("No index " + index +
                    " in Sequence instance of length " + size());
        AlignmentColumn column = get(index);
        return (SequenceAlignmentCell) column.cell(0);
    }

    public char getChar(int index) {
        return cell(index).getChar();

    }

    public String getCharAsString(int index) {
        return cell(index).getCharAsString();

    }

    public SecondaryStructure getSs(int index) {
        SecondaryStructure ss = (SecondaryStructure) cell(index).getAttribute(MeshiAttribute.SECONDARY_STRUCTURE_ATTRIBUTE);
        if (ss == null) return SecondaryStructure.secondaryStructure('C');
        return ss;
    }

    public char charAt(int index) {
        return cell(index).getChar();
    }

    public MeshiSequence renumber(MeshiSequence reference, AlignmentScheme alSchem) {//Tommer 23.9.14
        SequenceAlignment alignment = SequenceAlignment.substitutionAlignment(reference, this, alSchem, null);  // Not sure what this method is about. Chen 25.8.2017
        //May need updating to new constructor! Tommer 2.9.14
        SequenceAlignmentColumn firstColumn = alignment.get(0);
        if (reference.size() == 0) return null;
        if (!firstColumn.isExactMachWithGaps()) return renumber(reference.tail(), alSchem);
        if (!alignment.isExactMachWithGaps())
            throw new RuntimeException("Cannot renumber \n" + this + "\n" +
                    "with\n" +
                    alignment);
        return renumber(alignment);
    }

    public MeshiSequence tail() {
        MeshiSequence out = new MeshiSequence(comment(), charFilter);
        Iterator columns = iterator();
        if (!columns.hasNext()) return out;
        columns.next();
        while (columns.hasNext())
            out.add((SequenceAlignmentColumn) columns.next());
        return out;
    }

    public MeshiSequence renumber(SequenceAlignment alignment) {
        AlignmentCell pivot = null;
        int shift = 999999;
        for (int i = 0; (i < size()) & (pivot == null); i++) {
            SequenceAlignmentCell cell = cell(i);
            for (Iterator columns = alignment.iterator(); (columns.hasNext()) & (pivot == null);) {
                SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
                if ((column.cell(1) == cell) &
                        (!cell.gap())) {
                    pivot = column.cell(0);
                    shift = pivot.number - cell.number;
                }
            }
        }
        if (pivot == null)
            throw new RuntimeException("Cannot renumber \n" + this + "\n" +
                    "with\n" +
                    alignment);

        MeshiSequence out = new MeshiSequence(comment(), charFilter);
        for (int i = 0; (i < size()); i++) {
            SequenceAlignmentCell cell = cell(i);
            if (cell.gap()) out.add(cell);
            else out.add(new SequenceAlignmentCell(cell.getChar(), cell.number + shift));
        }
        return out;
    }

    public boolean add(SequenceAlignmentCell cell) {
        return add(new SequenceAlignmentColumn(cell));
    }

    private static class kolDichfin extends SequenceCharFilter {
        public boolean accept(Object obj) {
            return true;
        }
    }

    private static String getChainSequence(Chain chain) {
        String[] all = chain.toString().split("\n");
        return all[1];
    }

}
