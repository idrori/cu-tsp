/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.parameters.*;

/**
 * <pre>
 *              O
 *              |
 *     N - CA - C...n
 *     |
 *     H
 */
public class Gly extends ResidueExtendedAtoms {
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "             O\n" +
            "             |\n" +
            "    N - CA - C...n\n" +
            "    |\n" +
            "    H\n";

    public Gly(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.GLY, atomList, id, mode,molecularSystem);
    }

    public String comment() {
        return COMMENT;
    }
}

