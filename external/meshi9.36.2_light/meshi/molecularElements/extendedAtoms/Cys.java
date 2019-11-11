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
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB
 *           |
 *           SG
 */
public class Cys extends ResidueExtendedAtoms {
    public final Atom SG;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB\n" +
            "           |\n" +
            "          SG \n";

    public Cys(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.CYS, atomList, id, mode,molecularSystem);
        Object[] temp = new Object[1];
        SG = getAtom("SG", AtomType.CSG, atomList, this,molecularSystem);
        getAtoms().add(SG);
        bonds.add(CB.bond(SG));
    }

    public String comment() {
        return COMMENT;
    }
}
