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
 * ARGININ
 * <p/>
 *              O
 *              |
 *     N - CA - C...n
 *     |   |
 *     H   CB
 *         |
 *         CG - CD - NE - CZ - NH1
 *                        |
 *                        NH2
 */
public class Arg extends ResidueExtendedAtoms {
    public final Atom CG, CD, NE, HE, CZ, NH1, NH2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "              O\n" +
            "              |\n" +
            "     N - CA - C...n\n" +
            "     |   |\n" +
            "     H   CB\n" +
            "         |\n" +
            "         CG - CD - NE - CZ - NH1\n" +
            "                   |    |\n" +
            "                   HE   NH2\n";


    public Arg(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.ARG, atomList, id, mode, molecularSystem);

        getAtoms().add(CG = getAtom("CG", AtomType.RCG, atomList, this,molecularSystem));
        getAtoms().add(CD = getAtom("CD", AtomType.RCD, atomList, this,molecularSystem));
        getAtoms().add(NE = getAtom("NE", AtomType.RNE, atomList, this,molecularSystem));
        getAtoms().add(HE = getAtom("HE", AtomType.RHE, atomList, this,molecularSystem));
        getAtoms().add(CZ = getAtom("CZ", AtomType.RCZ, atomList, this,molecularSystem));
        getAtoms().add(NH1 = getAtom("NH1", AtomType.RNH, atomList, this,molecularSystem));
        getAtoms().add(NH2 = getAtom("NH2", AtomType.RNH, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(CD));
        bonds.add(CD.bond(NE));
        bonds.add(NE.bond(CZ));
        bonds.add(NE.bond(HE));
        bonds.add(CZ.bond(NH1));
        bonds.add(CZ.bond(NH2));
    }

    public void addHydrogens(BondParametersList bondParametersList, AngleParametersList angleParametersList) {
            super.addHydrogens(bondParametersList, angleParametersList);
            if (HE.nowhere() &&
                    (!CD.nowhere()) &&
                    (!NE.nowhere()) &&
                    (!CZ.nowhere())) {
                PutHpos.pos(HE, bondParametersList, angleParametersList);
            }
            //	    if ((NE != null) && HE.nowhere()) PutHpos.pos(HE,bondParametersList, angleParametersList);
    }

    public String comment() {
        return COMMENT;
    }
}
