package meshi.sequences;

import meshi.parameters.LocalStructure;

/**
 * Created by user on 15/02/2017.
 */
public class LocalStructureSequence extends ResidueSequence {



    public LocalStructureSequence(String sequence, LocalStructure[] localSsOfResidues, String comment) {
        super();
        Character weirdChar = weirdChar(sequence);
        if (sequence.length() != localSsOfResidues.length) throw new RuntimeException("Residue sequence and local structure length doesn't match. For more information: "+this.getClass().getName());
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
                cell.addAttribute(localSsOfResidues[i]);
                SequenceAlignmentColumn column = new SequenceAlignmentColumn(cell);
                add(column);
            }
        } else {
            throw new ResidueTypeException(weirdChar, sequence);
        }
    }

    private static class LocalStructureSequenceCharFilter extends SequenceCharFilter {
        public boolean accept(Object obj) {
            return true;
        }
    }
}
