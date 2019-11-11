/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;
import meshi.parameters.ResidueMode;
import meshi.parameters.ResidueType;

/**
 * <pre>
 * <p/>
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB
 *           |
 *           CG - OD2
 *           | (-)
 *           OD1
 */
public class Asp extends ResidueExtendedAtoms {
    public final Atom CG, OD1, OD2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB\n" +
            "           |\n" +
            "           CG - OD2\n" +
            "           | \n" +
            "           OD1\n";

    public Asp(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.ASP, atomList, id, mode,molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.DCG, atomList, this,molecularSystem));
        getAtoms().add(OD1 = getAtom("OD1", AtomType.DOD, atomList, this,molecularSystem));
        getAtoms().add(OD2 = getAtom("OD2", AtomType.DOD, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(OD2));
        bonds.add(CG.bond(OD1));
    }

    public String comment() {
        return COMMENT;
    }
}

