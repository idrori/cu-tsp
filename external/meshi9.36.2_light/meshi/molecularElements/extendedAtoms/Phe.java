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
 *           CG - CD1 - CE1
 *           |          |
 *           CD2- CE2 - CZ
 */
public class Phe extends ResidueExtendedAtoms {
    public final Atom CG, CD1, CE1, CZ, CD2, CE2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "               O\n" +
            "               |\n" +
            "      N - CA - C...n\n" +
            "      |   |\n" +
            "      H   CB\n" +
            "          |\n" +
            "          CG - CD1 - CE1\n" +
            "          |          |\n" +
            "          CD2- CE2 - CZ\n";


    public Phe(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.PHE, atomList, id, mode,molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.FCG, atomList, this,molecularSystem));
        getAtoms().add(CD1 = getAtom("CD1", AtomType.FCD, atomList, this,molecularSystem));
        getAtoms().add(CE1 = getAtom("CE1", AtomType.FCE, atomList, this,molecularSystem));
        getAtoms().add(CZ = getAtom("CZ", AtomType.FCZ, atomList, this,molecularSystem));
        getAtoms().add(CD2 = getAtom("CD2", AtomType.FCD, atomList, this,molecularSystem));
        getAtoms().add(CE2 = getAtom("CE2", AtomType.FCE, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(CD1));
        bonds.add(CG.bond(CD2));
        bonds.add(CD1.bond(CE1));
        bonds.add(CD2.bond(CE2));
        bonds.add(CE1.bond(CZ));
        bonds.add(CE2.bond(CZ));
    }

    public String comment() {
        return COMMENT;
    }
}
