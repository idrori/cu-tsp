/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.parameters;

import meshi.util.*;

public enum SecondaryStructure implements MeshiAttribute {
    HELIX("H"),
    SHEET("E"),
    COIL("C"),
    HELIX_OR_COIL("h"),
    SHEET_OR_COIL("e"),
    ALL("A"),
    UNK("X"),
    THREE10HELIX("G"),
    PIHELIX("S"),
    TURN("T"),
    GAP("-");

    public String getNameOneLetter() {
        return nameOneLetter;
    }

    private final String nameOneLetter;
    private static  double prob;

    private SecondaryStructure(String nameOneLetter) {
        this.nameOneLetter = nameOneLetter;
    }

    public static SecondaryStructure secondaryStructure(char c) {
    	if(c == '.')	c='-';
    	if((c >= 'a') && (c<='z'))	prob=0.5;
    	else if ((c >= 'A') && (c<='Z')) prob=1;
    	char newC = Character.toUpperCase(c);
        for (SecondaryStructure sec : SecondaryStructure.values())
            if (sec.nameOneLetter.charAt(0) == newC) return sec;
        throw new RuntimeException("Undefined secondary structure: " + newC + "\n" + "Please take a look at meshi.parameters.SecondaryStructure.");
    }

    public static boolean isSecondaryStructure(char c) {
        for (SecondaryStructure sec : SecondaryStructure.values())
            if ((sec.nameOneLetter.charAt(0) == c) &
                    (sec != ALL) &
                    (sec != UNK)) return true;
        return false;
    }

    public boolean equalsIgnorOr(SecondaryStructure other) {
        return ((this == HELIX && other == HELIX_OR_COIL) || (this == HELIX_OR_COIL && other == HELIX) ||
                (this == SHEET && other == SHEET_OR_COIL) || (this == SHEET_OR_COIL && other == SHEET) ||
                (this == COIL && other == HELIX_OR_COIL) || (this == HELIX_OR_COIL && other == COIL) ||
                (this == COIL && other == SHEET_OR_COIL) || (this == SHEET_OR_COIL && other == COIL));

    }

    public int key() {
        return SECONDARY_STRUCTURE_ATTRIBUTE;
    }
}
