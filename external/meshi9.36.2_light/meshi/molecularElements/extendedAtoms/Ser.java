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
 *                O
 *                |
 *       N - CA - C...n
 *           |
 *           CB
 *          /
 *        OG
 */
public class Ser extends ResidueExtendedAtoms {
    public final Atom OG;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "           |\n" +
            "           CB\n" +
            "          /\n" +
            "        OG\n";

    public Ser(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.SER, atomList, id, mode,molecularSystem);
        getAtoms().add(OG = getAtom("OG", AtomType.SOG, atomList, this,molecularSystem));
        bonds.add(CB.bond(OG));
    }

    public String comment() {
        return COMMENT;
    }
}

