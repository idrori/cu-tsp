/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction;

import meshi.util.*;
import meshi.util.filters.*;
import meshi.sequences.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

public class OriginalAtom implements MeshiAttribute {
    public static final Filter filter = new OriginalAtomFilter();
    public static final OriginalAtom attribute = new OriginalAtom();

    public int key() {
        return ORIGINAL_ATOM;
    }

    private static class OriginalAtomFilter implements Filter {
        public boolean accept(Object obj) {
            if (obj instanceof ResidueAlignmentColumn) {
                ResidueAlignmentColumn column = (ResidueAlignmentColumn) obj;
                Residue residue0 = column.residue0();
                Residue residue1 = column.residue1();
                if (residue0.dummy() | residue1.dummy()) return false;

                Atom ca0 = residue0.ca();
                Atom ca1 = residue1.ca();
                if (ca0.nowhere() | ca1.nowhere()) return false;

                OriginalAtom oa0 = (OriginalAtom) column.residue0().ca().getAttribute(ORIGINAL_ATOM);
                OriginalAtom oa1 = (OriginalAtom) column.residue1().ca().getAttribute(ORIGINAL_ATOM);
                if ((oa0 != null) | (oa1 != null)) return true;
                return false;
            } else {
                if (obj instanceof Atom) {
                    Atom atom = (Atom) obj;
                    OriginalAtom oa = (OriginalAtom) atom.getAttribute(ORIGINAL_ATOM);
                    return (oa != null);
                } else throw new RuntimeException("This filter operates only isOn getAtoms or ResidueAlignmentColumn");
            }
        }
    }

}