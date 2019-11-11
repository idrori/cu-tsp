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
 *           CG - SD - CE
 */
public class Met extends ResidueExtendedAtoms {
    public final Atom CG, SD, CE;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                 O\n" +
            "                 |\n" +
            "        N - CA - C...n\n" +
            "        |   |\n" +
            "        H   CB\n" +
            "            |\n" +
            "            CG - SD - CE\n";


    public Met(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.MET, atomList, id, mode,molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.MCG, atomList, this,molecularSystem));
        getAtoms().add(SD = getAtom("SD", AtomType.MSD, atomList, this,molecularSystem));
        getAtoms().add(CE = getAtom("CE", AtomType.MCE, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(SD));
        bonds.add(SD.bond(CE));
    }

    public String comment() {
        return COMMENT;
    }
}
