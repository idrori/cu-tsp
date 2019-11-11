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
 *       |   |
 *       CD  CB
 *       \  /
 *        CG
 */
public class Pro extends ResidueExtendedAtoms {
    public final Atom CG, CD;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       CD  CB\n" +
            "       \\  /\n" +
            "        CG\n";


    public Pro(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.PRO, atomList, id, mode,molecularSystem);
        getAtoms().add(CG = getAtom("CG", AtomType.PCG, atomList, this, molecularSystem));
        getAtoms().add(CD = getAtom("CD", AtomType.PCD, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(CD));
        bonds.add(CD.bond(N));
    }

    public String comment() {
        return COMMENT;
    }
}
