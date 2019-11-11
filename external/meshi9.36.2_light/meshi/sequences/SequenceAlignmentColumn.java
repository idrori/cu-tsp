/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import java.util.*;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.*;

public class SequenceAlignmentColumn extends AlignmentColumn implements MeshiAttribute {

    public SequenceAlignmentColumn(ResidueAlignmentColumn residueAlignmentColumn) {
        super(residueAlignmentColumn.cells.length);
        for (int i = 0; i < residueAlignmentColumn.cells.length; i++) {
            ResidueAlignmentCell oldCell = (ResidueAlignmentCell) residueAlignmentColumn.cell(i);
            Residue residue = oldCell.residue();
            String name = residue.type().nameOneLetter();
            char name1 = name.charAt(0);
            int number = residue.number();
            SequenceAlignmentCell newCell = new SequenceAlignmentCell(name1, number);
            add(i, newCell);
        }
    }

    public SequenceAlignmentColumn(int numberOfRows) {
        super(numberOfRows);
    }

    public SequenceAlignmentColumn(SequenceAlignmentCell cell) {
        super(1);
        add(0, cell);
    }

    public SequenceAlignmentColumn(SequenceAlignmentCell cell0, SequenceAlignmentCell cell1) {
        super(2);
        add(0, cell0);
        add(1, cell1);
    }

    public SequenceAlignmentColumn(SequenceAlignmentColumn sourceColumn) {
        super(sourceColumn.size());
        int i = 0;
        for (AlignmentCell sourceCell : sourceColumn.cells) {
            add(i, new SequenceAlignmentCell((SequenceAlignmentCell) sourceCell));
            i++;
        }
    }

    public void add(int index, char char1, int position) {
        add(index, createCell(char1, position));
    }

    public SequenceAlignmentColumn(char char0, int number0, char char1, int number1) {
        super(createCell(char0, number0), createCell(char1, number1));
        if ((char0 == char1) & (char0 == SequenceAlignmentCell.GAP_CHAR))
            throw new RuntimeException("Aligning two gaps is weird");
    }

    private static SequenceAlignmentCell createCell(char char1, int number1) {
        if (char1 == SequenceAlignmentCell.GAP_CHAR)
            return new SequenceAlignmentCell(number1);
        else return new SequenceAlignmentCell(char1, number1);
    }

    public char getChar(int index) {
        return ((SequenceAlignmentCell) cell(index)).getChar();
    }

    public boolean isExactMachWithGaps() {
        SequenceAlignmentCell cell;
        int nonGap = -1;
        for (int i = 0; i < cells.length; i++) {
            cell = (SequenceAlignmentCell) cell(i);
            if ((!cell.gap()) & (!cell.wildcard())) nonGap = i;
        }
        if (nonGap == -1) return true;
        for (int i = 0; i < cells.length; i++) {
            cell = (SequenceAlignmentCell) cell(i);
            if ((!cell.gap()) & (!cell.wildcard()) &
                    (!cell.equals(cell(nonGap)))) return false;
        }
        return true;
    }

    public boolean isExactMach() {
        int nonGap = -1;
        for (int i = 0; i < cells.length; i++)
            if (!cell(i).equals(cell(0))) return false;
        return true;
    }

    public int key() {
        return SEQUENCE_ALIGNMENT_COLUMN_ATTRIBUTE;
    }
    
    public String printLetter() {//Tommer 21.9.14
        String out = "";
        if (1 < cells.length) return "the method printLetter is to be used only with single cell S.A. Columns";
        out+=cells[0].obj;
        return out;
    }
}
		
	
