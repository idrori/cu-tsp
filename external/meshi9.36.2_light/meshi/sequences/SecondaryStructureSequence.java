/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import java.util.*;

public class SecondaryStructureSequence extends MeshiSequence {
    public static final SecondaryStructureSequenceCharFilter secondaryStructureCharFilter = new SecondaryStructureSequenceCharFilter();

    public SecondaryStructureSequence(String sequence, String comment) {
        super(sequence, comment, secondaryStructureCharFilter);
    }

    private static class SecondaryStructureSequenceCharFilter extends SequenceCharFilter {
        public boolean accept(Object obj) {
            Character c = ((Character) obj).charValue();
            if ("HEC".indexOf(c) >= 0) return true;
            return false;
        }
    }
}
