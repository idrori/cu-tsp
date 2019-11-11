/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

/**
 * A container for an ordered setResidue of coresponding protein elements.
 * Each of the elements comes form a different protein.
 */
public class AlignmentColumn {
    protected AlignmentCell[] cells;


    /**
     * Utility constructor for the sub-classes.
     */
    protected AlignmentColumn(int numberOfRows) {
        cells = new AlignmentCell[numberOfRows];
    }

    public AlignmentColumn(AlignmentCell cell0, AlignmentCell cell1) {
        cells = new AlignmentCell[2];
        cells[0] = cell0;
        cells[1] = cell1;
    }

    public void add(int index, AlignmentCell cell) {
        cells[index] = cell;
    }

    public AlignmentCell cell0() {
        return cell(0);
    }

    public AlignmentCell cell1() {
        return cell(1);
    }

    public AlignmentCell cell(int index) {
        if (cells[index] == null) throw new RuntimeException("empty cell");
        return cells[index];
    }

    public boolean hasGap() {
        int length = cells.length;
        for (int i = 0; i < length; i++)
            if (cells[i].gap()) return true;
        return false;
    }

    public boolean allGaps() {
        int length = cells.length;
        for (int i = 0; i < length; i++) {
            if (!cells[i].gap()) return false;
        }
        return true;
    }

    public String toString() {
        String out = "";
        for (int i = 0; i < cells.length; i++) {
            AlignmentCell cell = cells[i];
            out += cell.obj + "_" + cell.number + "\t";
        }
        return out;
    }

    public int size() {
        return cells.length;
    }
}
