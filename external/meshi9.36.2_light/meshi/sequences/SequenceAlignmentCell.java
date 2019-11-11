/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.util.*;

public class SequenceAlignmentCell extends AlignmentCell implements MeshiAttribute {
    public static final char GAP_CHAR = '-';
    public static final Character GAP = new Character(GAP_CHAR);
    public static final char WILDCARD_CHAR = 'X';
    public static final Character WILDCARD = new Character(WILDCARD_CHAR);

    public SequenceAlignmentCell(char c, int number, String comment) {

        super(new Character(c), number,comment);
    }
    public SequenceAlignmentCell(char c, int number) {
        super(new Character(c), number);
    }

    public SequenceAlignmentCell(int number) {
        super(GAP, number);
    }

    public SequenceAlignmentCell() {
        super(GAP, -1);
    }

    public SequenceAlignmentCell(SequenceAlignmentCell sourceCell) {
        super(new Character(sourceCell.getChar()), sourceCell.number());
    }

    public char getChar() {
        return ((Character) obj).charValue();
    }

    public String getCharAsString() {
        char[] cc = new char[1];
        cc[0] = getChar();
        return new String(cc);
    }


    public boolean gap() {
        return obj.equals(GAP);
    }

    public boolean wildcard() {
        return obj.equals(WILDCARD);
    }

    public boolean equals(Object obj) {
        SequenceAlignmentCell other = (SequenceAlignmentCell) obj;
        return (getChar() == other.getChar());
    }

    public static boolean wildcardOrGap(SequenceAlignmentCell cell) {
        if (cell.obj.equals(WILDCARD)) return true;
        if (cell.obj.equals(GAP)) return true;
        return false;
    }

    public static boolean wildcardOrGap(SequenceAlignmentColumn column) {
        for (int i = 0; i < column.size(); i++) {
            if (wildcardOrGap((SequenceAlignmentCell) column.cell(i))) return true;
        }
        return false;
    }

    public String toString() {
        return "" + getChar() + " " + number;
    }

    public MeshiAttribute getAttribute() {
        return getAttribute(key());
    }

    public int key() {
        return SEQUENCE_ALIGNMENT_COLUMN_ATTRIBUTE;
    }
}
    
