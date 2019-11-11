/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.geometry.putH.*;
import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.angle.*;

/**
 * <pre>
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB   CE3 - CZ3 - CH2
 *           |    |           |
 *           CG - CD2 - CE2 - CZ2
 *           |           |
 *           CD1 ------ NE1 - HE1
 */
public class Trp extends ResidueExtendedAtoms {
    public final Atom CG, CD1, NE1, HE1, CE2, CZ2, CH2, CZ3, CE3, CD2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB   CE3 - CZ3 - CH2\n" +
            "           |    |           |\n" +
            "           CG - CD2 - CE2 - CZ2\n" +
            "           |           |\n" +
            "           CD1 ------ NE1 - HE1\n";


    public Trp(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.TRP, atomList, id, mode,molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.WCG, atomList, this,molecularSystem));
        getAtoms().add(CD1 = getAtom("CD1", AtomType.WCD1, atomList, this,molecularSystem));
        getAtoms().add(NE1 = getAtom("NE1", AtomType.WNE, atomList, this,molecularSystem));
        getAtoms().add(HE1 = getAtom("HE1", AtomType.WHE, atomList, this,molecularSystem));
        getAtoms().add(CD2 = getAtom("CD2", AtomType.WCD2, atomList, this,molecularSystem));
        getAtoms().add(CE2 = getAtom("CE2", AtomType.WCE2, atomList, this,molecularSystem));
        getAtoms().add(CE3 = getAtom("CE3", AtomType.WCE3, atomList, this,molecularSystem));

        getAtoms().add(CZ3 = getAtom("CZ3", AtomType.WCZ3, atomList, this,molecularSystem));
        getAtoms().add(CH2 = getAtom("CH2", AtomType.WCH2, atomList, this,molecularSystem));
        getAtoms().add(CZ2 = getAtom("CZ2", AtomType.WCZ2, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(CD1));
        bonds.add(CG.bond(CD2));
        bonds.add(CD1.bond(NE1));
        bonds.add(CD2.bond(CE2));
        bonds.add(NE1.bond(HE1));
        bonds.add(CE2.bond(CZ2));
        bonds.add(CZ2.bond(CH2));
        bonds.add(CH2.bond(CZ3));
        bonds.add(CZ3.bond(CE3));
        bonds.add(CE3.bond(CD2));
        bonds.add(NE1.bond(CE2));
    }

    public String comment() {
        return COMMENT;
    }

    public void addHydrogens(BondParametersList bondParametersList, AngleParametersList angleParametersList) {
        super.addHydrogens(bondParametersList, angleParametersList);
            if (HE1.nowhere() &&
                    (!CD1.nowhere()) &&
                    (!NE1.nowhere()) &&
                    (!CE2.nowhere())) PutHpos.pos(HE1, bondParametersList, angleParametersList);
    }
}
