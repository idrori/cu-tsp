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
 *           |
 *           CB
 *          /  \
 *        OG1  CG2
 */
public class Thr extends ResidueExtendedAtoms {
    public final Atom OG1, CG2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "           |\n" +
            "           CB\n" +
            "          /  \\\n" +
            "        OG1  CG2\n";


    public Thr(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.THR, atomList, id, mode,molecularSystem);

        getAtoms().add(CG2 = getAtom("CG2", AtomType.TCG, atomList, this,molecularSystem));
        getAtoms().add(OG1 = getAtom("OG1", AtomType.TOG, atomList, this,molecularSystem));
        bonds.add(CB.bond(OG1));
        bonds.add(CB.bond(CG2));
    }

    public String comment() {
        return COMMENT;
    }
}
