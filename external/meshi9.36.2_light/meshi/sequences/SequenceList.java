/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import java.util.*;

import meshi.util.*;
import meshi.util.file.*;
import meshi.util.filters.Filter;

public class SequenceList extends ArrayList<MeshiSequence> implements KeyWords {
    private String fileName = null;
    protected static final SequenceCharFilter acceptAllChars = new AcceptAll();

    public SequenceList() {
        super();
    }

    public SequenceList(String fileName) {
        if (!FastaList.isFasta(fileName)) throw new RuntimeException("currently only Fasta files are used");
        SequenceList temp = new FastaList(fileName);
        for (MeshiSequence seq : temp)
            add(seq);
        this.fileName = fileName;
    }

    public SequenceList(SequenceAlignment alignment) {
        super();
        SequenceAlignmentColumn firstColumn = alignment.get(0);
        int numberOfRows = firstColumn.size();
        for (int i = 0; i < numberOfRows; i++) {
            MeshiSequence sequence = new MeshiSequence(alignment.comments.get(i),
                    acceptAllChars);
            for (Iterator columns = alignment.iterator(); columns.hasNext();) {
                SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
                sequence.add((SequenceAlignmentCell) column.cell(i));
            }
            add(sequence);
        }
    }

    public SequenceList(SequenceList source) {
        super();
        for (MeshiSequence sourceSequence : source) {
            add(new MeshiSequence(sourceSequence));
        }
    }


    //------------------------ extract sequences from the list ----------------------------
    private MeshiSequence getSequence(String name) {
        for (Iterator sequences = iterator(); sequences.hasNext();) {
            MeshiSequence sequence = (MeshiSequence) sequences.next();
            if (sequence.comment().indexOf(name) > -1) return sequence;
        }
        return null;
    }

    public ResidueSequence getResidueSequence(String name) {
        MeshiSequence sequence = getSequence(name);
        if (sequence == null) return null;
        return new ResidueSequence(sequence);
    }

    public ResidueSequence getResidueSequence() {
        return getResidueSequence(AA_SEQUENCE.key);
    }

    public ResidueSequence getResidueSequence(CommandList commands, Key seqNameKey) {
        String name = commands.firstWord(seqNameKey).secondWord();
        return getResidueSequence(name);
    }


    private static class IsSequenceFilter implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof MeshiSequence);
        }
    }

    public String toString() {
        String out = "Sequence list \n";
        for (Iterator sequences = iterator(); sequences.hasNext();) {
            out += sequences.next().toString() + "\n";
        }
        return out;
    }

    public String fileName() {
        return fileName;
    }

    private static class AcceptAll extends SequenceCharFilter {
        public boolean accept(Object obj) {
            Character c = ((Character) obj).charValue();
            return true;
        }

    }

    public void print() {
        for (MeshiSequence s : this)
            System.out.println(s);
    }

}
	
