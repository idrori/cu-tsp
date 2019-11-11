/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import java.util.*;

public class AccesibilitySequence extends MeshiSequence {
    public static final AccesibilitySequenceCharFilter AccesibilityCharFilter = new AccesibilitySequenceCharFilter();

    public AccesibilitySequence(String[] accesibility, String comment) {
        super(accesibility, comment, AccesibilityCharFilter);
    }

    private static class AccesibilitySequenceCharFilter extends SequenceCharFilter {
        public boolean accept(Object obj) {
            Character c = ((Character) obj).charValue();
            if ("AB".indexOf(c) >= 0) return true;
            return false;
        }
    }
}
