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
 *       H   CB
 *           |
 *           CG - CD - NE2 - HE21
 *                |    |
 *                OE1  HE22
 */
public class Gln extends ResidueExtendedAtoms {
    public final Atom CG, CD, OE1, NE2, HE21, HE22;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB\n" +
            "           |\n" +
            "           CG - CD - NE2 - HE21\n" +
            "                |    |\n" +
            "                OE1  HE22\n";

    public Gln(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.GLN, atomList, id, mode,molecularSystem);
        getAtoms().add(CG = getAtom("CG", AtomType.QCG, atomList, this,molecularSystem));
        getAtoms().add(CD = getAtom("CD", AtomType.QCD, atomList, this,molecularSystem));
        getAtoms().add(NE2 = getAtom("NE2", AtomType.QNE, atomList, this,molecularSystem));
        getAtoms().add(HE21 = getAtom("HE21", AtomType.QHE1, atomList, this,molecularSystem));
        getAtoms().add(HE22 = getAtom("HE22", AtomType.QHE2, atomList, this, molecularSystem));
        getAtoms().add(OE1 = getAtom("OE1", AtomType.QOE, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(CD));
        bonds.add(CD.bond(OE1));
        bonds.add(CD.bond(NE2));
        bonds.add(NE2.bond(HE21));
        bonds.add(NE2.bond(HE22));
    }

    public String comment() {
        return COMMENT;
    }

    public void addHydrogens(BondParametersList bondParametersList, AngleParametersList angleParametersList) {
        super.addHydrogens(bondParametersList, angleParametersList);
            double minDis = 1000;
            Coordinates minCoor = new Coordinates(1000, 1000, 1000);
            if (HE21.nowhere() &&
                    (!NE2.nowhere()) &&
                    (!OE1.nowhere()) &&
                    (!CD.nowhere())) {
                if (!HE22.nowhere()) {
                    PutHpos.pos(HE21, bondParametersList, angleParametersList);
                } else {
                    for (int i = 0; i < 1000; i++) {
                        PutHpos.pos(HE21, bondParametersList, angleParametersList);
                        double dis = OE1.distanceFrom(HE21);
                        if (dis < minDis) {
                            minDis = dis;
                            minCoor.set(new Coordinates(HE21));
                        }
                        HE21.resetCoordinates();
                    }
                    HE21.setXYZ(minCoor.x(), minCoor.y(), minCoor.z());
                }
            }
            minDis = 1000;
            minCoor = new Coordinates(1000, 1000, 1000);
            if (HE22.nowhere() &&
                    (!NE2.nowhere()) &&
                    (!OE1.nowhere()) &&
                    (!CD.nowhere())) {
                if (!HE21.nowhere()) {
                    PutHpos.pos(HE22, bondParametersList, angleParametersList);
                } else {
                    for (int i = 0; i < 1000; i++) {
                        PutHpos.pos(HE22, bondParametersList, angleParametersList);
                        double dis = OE1.distanceFrom(HE22);
                        if (dis < minDis) {
                            minDis = dis;
                            minCoor.set(new Coordinates(HE22));
                        }
                        HE22.resetCoordinates();
                    }
                    HE22.setXYZ(minCoor.x(), minCoor.y(), minCoor.z());
                }
            }

    }
}
