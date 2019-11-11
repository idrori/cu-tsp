/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.sequences.aligner.*;
import meshi.util.string.StringList;

import java.util.ArrayList;
import java.util.Iterator;

public class SequenceAlignment extends ArrayList<SequenceAlignmentColumn> {
    private Double score = null;
    public final StringList comments;

    public SequenceAlignment() {
        super();
        comments = new StringList();
    }

    public SequenceAlignment(ResidueAlignment residueAlignment) {
        this();
        setScore(residueAlignment.score());
        for (Iterator resColumns = residueAlignment.iterator(); resColumns.hasNext();)
            add(new SequenceAlignmentColumn((ResidueAlignmentColumn) resColumns.next()));
        for (String s : residueAlignment.comments)
            comments.add(s);
    }

    public SequenceAlignment(MeshiSequence sequence1, MeshiSequence sequence2) {
        this();
        SequenceList sequenceList = new SequenceList();
        sequenceList.add(sequence1);
        sequenceList.add(sequence2);
        SequenceAlignment temp = new SequenceAlignment(sequenceList);
        for (SequenceAlignmentColumn sac : temp)
            add(sac);
        for (String s : temp.comments)
            comments.add(s);
    }

    public SequenceAlignment(SequenceList sequenceList) {
        this();
        int numberOfSequences = sequenceList.size();
        boolean done = false;

        int length = (sequenceList.get(0)).size();
        for (int i = 0; i < numberOfSequences; i++) {
            MeshiSequence sequence = sequenceList.get(i);
            if (sequence.size() != length) {
                System.out.println("A problem with sequences:");
                sequenceList.print();
                throw new RuntimeException("All sequences must have the same Length");
            }
            comments.add(sequence.comment());
        }

        for (int iPosition = 0; iPosition < length; iPosition++) {
            SequenceAlignmentColumn column = new SequenceAlignmentColumn(numberOfSequences);
            for (int iSequence = 0; iSequence < numberOfSequences; iSequence++) {
                MeshiSequence sequence = sequenceList.get(iSequence);
                column.add(iSequence, sequence.cell(iPosition));
            }
            add(column);
        }
    }


    public String toString() {
       // if (1 == 1) throw new RuntimeException("XXXXXXXXXXXXX") ;
        SequenceList sequenceList = new SequenceList(this);
        String out = "Score: "+this.score()+"\n";//Added by Tommer 7.9.14
        for (Iterator sequences = sequenceList.iterator(); sequences.hasNext();)
            out += "" + sequences.next();
        return out;
    }

    public boolean isExactMach() {
        for (Iterator columns = iterator(); columns.hasNext();)
            if (!((SequenceAlignmentColumn) columns.next()).isExactMach()) return false;
        return true;
    }

    public boolean isExactMachWithGaps() {
        for (Iterator columns = iterator(); columns.hasNext();)
            if (!((SequenceAlignmentColumn) columns.next()).isExactMachWithGaps()) return false;
        return true;
    }

    public static SequenceAlignment identityAlignment(MeshiSequence sequence1, MeshiSequence sequence2) {
        return substitutionAlignment(sequence1,sequence2,new IdentityMatrix(), ResidueAlignmentMethod.IDENTITY);
    }

    public static SequenceAlignment substitutionAlignment(MeshiSequence sequence1, MeshiSequence sequence2,
    		SubstitutionMatrix substitutionMatrix, ResidueAlignmentMethod method) {//name updated by Tommer 9.914, throws added 11.9.14
    	if(substitutionMatrix==null|| substitutionMatrix.order()==null)// exception condition need POLISHING Tommer 23.9.14
    		throw new RuntimeException("the code at this point needs to be adapted"
    				+ "to new identitiyAlignment method signature and provide a matrix");
    	//signature changed for identityAlignment Tommer 1.9.14
        DpMatrix matrix = new DpMatrix(sequence1, sequence2, new SubstitutionScorer(substitutionMatrix) ,substitutionMatrix.minScore());//Edited by Tommer 1.9.14, gap penalty was -0.2
      //signature changed for Identity Tommer 1.9.14
        SequenceAlignment out = matrix.backTrack(method);

        return out;
    }

    public double score() {
        if (score == null) throw new RuntimeException("Alignment without score");
        return score.doubleValue();
    }

    public void setScore(double score) {
        this.score = new Double(score);
    }

    public void print() {
        for (SequenceAlignmentColumn sac : this)
            System.out.println(sac);
    }
    
    public void printAsLine() {//Tommer 21.9.14
        for (SequenceAlignmentColumn sac : this)
            System.out.print(sac.printLetter());
        System.out.println("");
    }

    /**
     * Fetches a column with a specific cell number in a specific row. This column now is
     * a handle to the corresponding elements of the other proteins.
     */
    public AlignmentColumn getColumn(int row, int number) {
        for (Iterator columns = iterator(); columns.hasNext();) {
            AlignmentColumn column = (AlignmentColumn) columns.next();
            if (column.cell(row).number == number) return column;
        }
        return null;
    }

    public int countIdentities() {
        int out = 0;
        for (AlignmentColumn column : this) {
            AlignmentCell cell0 = column.cell0();
            AlignmentCell cell1 = column.cell1();
            if (cell0.object().equals(cell1.object()))
                out++;
        }
        return out;
    }

    public void addAll(SequenceAlignment sequenceAlignment) {
        super.addAll(sequenceAlignment);
        for (String s : sequenceAlignment.comments)
            comments.add(s);
        score = sequenceAlignment.score();
    }
}
