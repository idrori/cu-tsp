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
 *       H   CB
 *           |
 *           CG - CD - CE - NZ
 */
public class Lys extends ResidueExtendedAtoms {
    public final Atom CG, CD, CE, NZ;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB\n" +
            "           |\n" +
            "           CG - CD - CE - NZ\n";

    public Lys(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.LYS, atomList, id, mode,molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.KCG, atomList, this,molecularSystem));
        getAtoms().add(CD = getAtom("CD", AtomType.KCD, atomList, this,molecularSystem));
        getAtoms().add(CE = getAtom("CE", AtomType.KCE, atomList, this,molecularSystem));
        getAtoms().add(NZ = getAtom("NZ", AtomType.KNZ, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(CD));
        bonds.add(CD.bond(CE));
        bonds.add(CE.bond(NZ));
    }

    public String comment() {
        return COMMENT;
    }
}
