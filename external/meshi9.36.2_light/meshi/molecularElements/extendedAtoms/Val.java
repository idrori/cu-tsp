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
 *       H   CB - CG1
 *           |
 *           CG2
 */
public class Val extends ResidueExtendedAtoms {
    public final Atom CG1, CG2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB - CG1\n" +
            "           |\n" +
            "           CG2\n";


    public Val(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.VAL, atomList, id, mode,molecularSystem);
        getAtoms().add(CG1 = getAtom("CG1", AtomType.VCG1, atomList, this,molecularSystem));
        getAtoms().add(CG2 = getAtom("CG2", AtomType.VCG2, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG1));
        bonds.add(CB.bond(CG2));
    }

    public String comment() {
        return COMMENT;
    }
}
