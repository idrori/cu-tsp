/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.util.*;
import meshi.util.file.*;
import meshi.util.filters.*;
import meshi.util.string.*;

import java.util.*;

public class FastaList extends SequenceList {
    public FastaList(String fileName) {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String sequence, comment, line;
        boolean empty;
        try {
            StringList lineList = new StringList(reader);
            for (Iterator lines = lineList.iterator(); lines.hasNext();) {
                comment = sequence = "";
                while (lines.hasNext() & comment.equals("")) {
                    line = (String) lines.next();
                    empty = (line.trim().equals("") | line.startsWith("#"));
                    if (!empty) {
                        if (line.startsWith(">"))
                            comment = line.substring(1).trim();
                        else throw new RuntimeException("unexpected line in " + fileName + "\n" +
                                line + "\n" +
                                "expected a line that start with >");
                    }
                }
                if (!comment.equals("")) {
                    while (lines.hasNext() & (!sequence.endsWith("*"))) {
                        line = (String) lines.next();
                        empty = (line.trim().equals("") | line.startsWith("#"));
                        if (!empty) {
                            sequence += removeSpaces(line.trim());
                        }
                    }
                    if (!sequence.endsWith("*")) {
                        throw new RuntimeException(" A Fasta sequence mast end in *");
                    }
                    add(new ResidueSequence(sequence.substring(0, sequence.length() - 1), // to remove the last character '*'.
                            comment));
                }
            }
        }
        catch (ResidueTypeException ex) {
            throw ex;
        }
        catch (Exception ex) {
            throw new RuntimeException("Bad Fasta file\n" + ex);
        }
    }

    public static boolean isFasta(String fileName) {
        try {
            new FastaList(fileName);
        }
        catch (ResidueTypeException ex) {
            throw new RuntimeException("\n" + "A problem in Fasta file - " + fileName + "\n" + ex);
        }
        catch (Exception ex) {
            System.out.println(ex);
            return false;
        }
        return true;
    }


    public static String removeSpaces(String source) {
        if (source.indexOf(' ') == -1) return source;
        String out = "";
        for (int index = 0; index < source.length(); index++) {
            String temp = source.substring(index, index - 1);
            if (!temp.equals(" ")) out += temp;
        }
        return out;
    }

}
	
