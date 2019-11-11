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
 *       H   CB - CG - CD1
 *                |
 *                CD2
 */
public class Leu extends ResidueExtendedAtoms {
    public final Atom CG, CD1, CD2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB - CG - CD1\n" +
            "                |\n" +
            "                CD2\n";


    public Leu(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.LEU, atomList, id, mode,molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.LCG, atomList, this,molecularSystem));
        getAtoms().add(CD1 = getAtom("CD1", AtomType.LCD1, atomList, this,molecularSystem));
        getAtoms().add(CD2 = getAtom("CD2", AtomType.LCD2, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(CD1));
        bonds.add(CG.bond(CD2));
    }

    public String comment() {
        return COMMENT;
    }
}
