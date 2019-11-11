/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.geometry.*;

/**
 * <pre>
 * Alanin from Levitt, JMB 168:592 (1983) table 2.
 *           CB  O
 *           |   |
 *      N - CA - C...n
 *      |
 *      H
 */
public class Ala extends ResidueExtendedAtoms {
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "      CB  O\n" +
            "      |   |\n" +
            "  N - CA - C...n\n" +
            "  |\n" +
            "  H\n";

    public Ala(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.ALA, atomList, id, mode,molecularSystem);
    }

    public String comment() {
        return COMMENT;
    }
}
