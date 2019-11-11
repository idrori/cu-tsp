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
 *                 O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB
 *           |
 *           CG - ND2 - HD21
 *           |    |
 *           OD1  HD22
 */
public class Asn extends ResidueExtendedAtoms {
    public final Atom CG, OD1, ND2, HD21, HD22;
    public static final String NAME = "ASN";
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n" +
            "                 O\n" +
            "                |\n" +
            "       N - CA - C...n\n" +
            "       |   |\n" +
            "       H   CB\n" +
            "           |\n" +
            "          CG - ND2 - HD21\n" +
            "           |    |\n" +
            "           OD1  HD22\n";


    public Asn(AtomList atomList, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        super(ResidueType.ASN, atomList, id, mode,molecularSystem);
        getAtoms().add(CG = getAtom("CG", AtomType.NCG, atomList, this,molecularSystem));
        getAtoms().add(ND2 = getAtom("ND2", AtomType.NND, atomList, this,molecularSystem));
        getAtoms().add(HD21 = getAtom("HD21", AtomType.NHD1, atomList, this,molecularSystem));
        getAtoms().add(HD22 = getAtom("HD22", AtomType.NHD2, atomList, this,molecularSystem));
        getAtoms().add(OD1 = getAtom("OD1", AtomType.NOD, atomList, this,molecularSystem));
        bonds.add(CB.bond(CG));
        bonds.add(CG.bond(OD1));
        bonds.add(CG.bond(ND2));
        bonds.add(ND2.bond(HD21));
        bonds.add(ND2.bond(HD22));
    }

    public String comment() {
        return COMMENT;
    }

    public void addHydrogens(BondParametersList bondParametersList, AngleParametersList angleParametersList) {
        super.addHydrogens(bondParametersList, angleParametersList);
            double minDis = 1000;
            Coordinates minCoor = new Coordinates(1000, 1000, 1000);
            if (HD21.nowhere() &&
                    (!ND2.nowhere()) &&
                    (!OD1.nowhere()) &&
                    (!CG.nowhere())) {
                if (!HD22.nowhere()) {
                    PutHpos.pos(HD21, bondParametersList, angleParametersList);
                } else {
                    for (int i = 0; i < 1000; i++) {
                        PutHpos.pos(HD21, bondParametersList, angleParametersList);
                        double dis = OD1.distanceFrom(HD21);
                        if (dis < minDis) {
                            minDis = dis;
                            minCoor.set(new Coordinates(HD21));
                        }
                        HD21.resetCoordinates();
                    }
                    HD21.setXYZ(minCoor.x(), minCoor.y(), minCoor.z());
                }
            }
            minDis = 1000;
            minCoor = new Coordinates(1000, 1000, 1000);
            if (HD22.nowhere() &&
                    (!ND2.nowhere()) &&
                    (!OD1.nowhere()) &&
                    (!CG.nowhere())) {
                if (!HD21.nowhere()) {
                    PutHpos.pos(HD22, bondParametersList, angleParametersList);
                } else {
                    for (int i = 0; i < 1000; i++) {
                        PutHpos.pos(HD22, bondParametersList, angleParametersList);
                        double dis = OD1.distanceFrom(HD22);
                        if (dis < minDis) {
                            minDis = dis;
                            minCoor.set(new Coordinates(HD22));
                        }
                        HD22.resetCoordinates();
                    }
                    HD22.setXYZ(minCoor.x(), minCoor.y(), minCoor.z());
                }
            }
    }
}
