/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.sequences.*;

public class OriginalAtoms implements Filter {
    Atom[] atomsArray;
    int size;

    public OriginalAtoms(AtomList atomList) {
        size = atomList.size();
        atomsArray = new Atom[size];
        for (int i = 0; i < size; i++) {
            atomsArray[i] = atomList.atomAt(i);
        }
    }

    public OriginalAtoms(ResidueAlignment residueAlignment, int row) {
        this(residueAlignment.getCaList(row));
    }

    public boolean accept(Object obj) {
        Atom atom = (Atom) obj;
        if (atom == null) throw new RuntimeException("(atom == null)");
        if (atomsArray == null) throw new RuntimeException("(atomsArray == null)");
        for (int i = 0; i < size; i++) {
            String name = atom.name();
            Atom iAtom = atomsArray[i];

            if (name.equals(iAtom.name()) &
                    (atom.residueNumber() == iAtom.residueNumber()) &
                    atom.residueName().equals(iAtom.residueName())) return true;
        }
        return false;
    }
}
