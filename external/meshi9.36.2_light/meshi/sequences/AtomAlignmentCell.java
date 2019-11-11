/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

public class AtomAlignmentCell extends AlignmentCell {
    public AtomAlignmentCell(Atom atom, String comment) {
        super(atom, atom.number(), comment);
    }

    public Atom atom() {
        return (Atom) obj;
    }

    public boolean gap() {
        return (obj == null);
    }

    public String toString() {
        return "{ " + atom() + " (" + atom().molecularSystem + ") }";
    }
}
