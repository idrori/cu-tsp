/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.geometry.putH.NotEnoughBoundAtomsException;
import meshi.geometry.putH.PutHpos;
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
 *           CG - CD2 - NE2 - HE2
 *           |           |
 *     HD1 - ND1 ------ CE1
 */
public class His extends ResidueExtendedAtoms {
    public final Atom CG, CD2, NE2, HE2, CE1, ND1, HD1;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "               O\n" +
            "               |\n" +
            "      N - CA - C...n\n" +
            "      |   |\n" +
            "      H   CB\n" +
            "          |\n" +
            "          CG - CD2 - NE2 - HE2\n" +
            "          |           |\n" +
            "    HD1 - ND1 ------ CE1\n";


    public His(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.HIS, atomList, id, mode,molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.HCG, atomList, this,molecularSystem));
        getAtoms().add(ND1 = getAtom("ND1", AtomType.HND, atomList, this,molecularSystem));
        getAtoms().add(HD1 = getAtom("HD1", AtomType.HHD, atomList, this,molecularSystem));
        getAtoms().add(CE1 = getAtom("CE1", AtomType.HCE, atomList, this,molecularSystem));
        getAtoms().add(NE2 = getAtom("NE2", AtomType.HNE, atomList, this,molecularSystem));
        getAtoms().add(HE2 = getAtom("HE2", AtomType.HHE, atomList, this,molecularSystem));
        getAtoms().add(CD2 = getAtom("CD2", AtomType.HCD, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(ND1));
        bonds.add(ND1.bond(HD1));
        bonds.add(ND1.bond(CE1));
        bonds.add(CG.bond(CD2));
        bonds.add(CD2.bond(NE2));
        bonds.add(CE1.bond(NE2));
        bonds.add(NE2.bond(HE2));
    }

    public String comment() {
        return COMMENT;
    }

    public void addHydrogens(BondParametersList bondParametersList, AngleParametersList angleParametersList) {
        super.addHydrogens(bondParametersList, angleParametersList);
            if (HD1.nowhere() &&
                    (!CG.nowhere()) &&
                    (!ND1.nowhere()) &&
                    (!CE1.nowhere())) {
                PutHpos.pos(HD1, bondParametersList, angleParametersList);
            }
            if (HE2.nowhere() &&
                    (!CD2.nowhere()) &&
                    (!NE2.nowhere()) &&
                    (!CE1.nowhere())) {
                PutHpos.pos(HE2, bondParametersList, angleParametersList);
            }
    }
}
