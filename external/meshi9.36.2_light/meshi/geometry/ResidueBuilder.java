/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

import java.util.*;


public class ResidueBuilder {

    private static double[] xyza = new double[4];

    public ResidueBuilder() {
    }

    public static void build(Residue res, int typeIndex, double[] chi1to4) {
        build(res, getType(typeIndex), chi1to4);
    }

    private static ResidueType getType(int index) {
        for (ResidueType type : ResidueType.values())
            if (type.ordinal() == index) return type;
        return null;
    }

    public static void build(Residue res, ResidueType type, double[] chi1to4) {
        if (res.dummy() || res.ca().nowhere()) return;
        switch (type) {

            case ALA:
                buildALA(res);//A
                break;
            case CYS:
                buildCYS(res, chi1to4);//C
                break;
            case ASP:
                buildASP(res, chi1to4);//D
                break;
            case GLU:
                buildGLU(res, chi1to4);//E
                break;
            case PHE:
                buildPHE(res, chi1to4);//F
                break;
            case GLY:
                buildGLY(res);//G
                break;
            case HIS:
                buildHIS(res, chi1to4);//H
                break;
            case ILE:
                buildILE(res, chi1to4);//I
                break;
            case LYS:
                buildLYS(res, chi1to4);//K
                break;
            case LEU:
                buildLEU(res, chi1to4);//L
                break;
            case MET:
                buildMET(res, chi1to4);//M
                break;
            case ASN:
                buildASN(res, chi1to4);//N
                break;
            case PRO:
                buildPRO(res, chi1to4);//P
                break;
            case GLN:
                buildGLN(res, chi1to4);//Q
                break;
            case ARG:
                buildARG(res, chi1to4);//R
                break;
            case SER:
                buildSER(res, chi1to4);//S
                break;
            case THR:
                buildTHR(res, chi1to4);//T
                break;
            case VAL:
                buildVAL(res, chi1to4);//V
                break;
            case TRP:
                buildTRP(res, chi1to4);//W
                break;
            case TYR:
                buildTYR(res, chi1to4);//Y
                break;
            default:
                throw new RuntimeException("wrong residue type: " + type);
        }
    }

    protected static void getAtom_xyza(double[] xyza, double bond, double angle, double chi, Atom A1, Atom A2, Atom A3) {
        if (A1 == null)
            throw new RuntimeException("Missing atom 1 in the torsion construction");
        if (A2 == null)
            throw new RuntimeException("Missing atom 2 in the torsion construction");
        if (A3 == null)
            throw new RuntimeException("Missing atom 3 in the torsion construction");
        /*xyza should be of type double[4]*/
        for (int i = 0; i < 4; i++)
            xyza[i] = 0;
        double moveTo[] = new double[4];
        moveTo[0] = bond * Math.sin(angle) * Math.sin(chi);
        moveTo[1] = bond * Math.sin(angle) * Math.cos(chi);
        moveTo[2] = bond * Math.cos(angle);
        moveTo[3] = 1;

        double M[][] = new double[4][4];
        double invM[][] = new double[4][4];
        ViewAt.transformToOrigin(M, invM, A1, A2, A3);
        // moveTo*invM
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                xyza[i] += moveTo[j] * invM[j][i];
            }
        }
    }


    public static void buildBackbone(Residue originalRes, double phi, double psi) {
        // Build O
        double bond = 1.233;
        double angle = (Math.PI / 180) * 120.54;

        getAtom_xyza(xyza, bond, angle, psi + Math.PI,
                getAtom(originalRes.getAtoms(), "C"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "O").setXYZ(xyza[0], xyza[1], xyza[2]);

        // Build H
        if (getAtom(originalRes.getAtoms(), "H") == null)
            return;
        bond = 1.0;
        angle = (Math.PI / 180) * 119.0;

        getAtom_xyza(xyza, bond, angle, phi + Math.PI,
                getAtom(originalRes.getAtoms(), "N"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "C"));

        getAtom(originalRes.getAtoms(), "H").setXYZ(xyza[0], xyza[1], xyza[2]);

    }

    protected static void buildALA(Residue originalRes) {
        // Build CB
        double bond = 1.54;
        double angle = (Math.PI / 180) * 109.8;

        getAtom_xyza(xyza, bond, angle, 122.2 * (Math.PI / 180),
                getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "C"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CB").setXYZ(xyza[0], xyza[1], xyza[2]);
    }


    protected static void buildCYS(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        //create SG
        double bond = 1.808;     /*        CH2E-SG1E */
        double angle = (Math.PI / 180) * 114.4; /*  CH1E-CH2E-SG1E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "SG").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildASP(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom ASP CG C */
        double bond = 1.516;     /*       CH2E-C */
        double angle = (Math.PI / 180) * 112.6; /*  CH1E-CH2E-C */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ASP OD1 OC */
        bond = 1.249;     /*       C-OC */
        angle = (Math.PI / 180) * 118.4;/*  CH2E-C-OC */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "OD1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ASP OD2 OC */
        bond = 1.249;     /*       C-OC */
        angle = (Math.PI / 180) * 118.4;/*  CH2E-C-OC */

        getAtom_xyza(xyza, bond, angle, chi[1] - Math.PI,
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "OD2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildGLU(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom GLU CG CH2E */
        double bond = 1.520;     /*       CH2E-CH2E */
        double angle = (Math.PI / 180) * 114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom GLU CD C */
        bond = 1.516;     /*       CH2E-C */
        angle = (Math.PI / 180) * 112.6; /*  CH2E-CH2E-C */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom GLU OE1 OC */
        bond = 1.249;     /*       C-OC */
        angle = (Math.PI / 180) * 118.4; /*  CH2E-C-OC */

        getAtom_xyza(xyza, bond, angle, chi[2],
                getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "OE1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom GLU OE2 OC */
        bond = 1.249;     /*       C-OC */
        angle = (Math.PI / 180) * 118.4; /*  CH2E-C-OC */

        getAtom_xyza(xyza, bond, angle, chi[2] + Math.PI,
                getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "OE2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildPHE(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom PHE CG CF */
        double bond = 1.502;     /*       CH2E-CF */
        double angle = (Math.PI / 180) * 113.8; /*  CH1E-CH2E-CF */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom PHE CD1 CR1E */
        bond = 1.384;     /*       CF-CR1E */
        angle = (Math.PI / 180) * 120.7; /*  CH2E-CF-CR1E */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom PHE CD2 CR1E */
        bond = 1.384;     /*       CF-CR1E */
        angle = (Math.PI / 180) * 120.7; /*  CH2E-CF-CR1E */

        getAtom_xyza(xyza, bond, angle, chi[1] - Math.PI,
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom PHE CE1 CR1E */
        bond = 1.382;     /*     CR1E-CR1E */
        angle = (Math.PI / 180) * 120.7; /*  CF-CR1E-CR1E */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CD1"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom PHE CE2 CR1E */
        bond = 1.382;     /*     CR1E-CR1E */
        angle = (Math.PI / 180) * 120.7; /*  CF-CR1E-CR1E */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom PHE CZ CR1E */
        bond = 1.382;     /*       CR1E-CR1E */
        angle = (Math.PI / 180) * 120.0; /*  CR1E-CR1E-CR1E */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "CE2"), getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "CZ").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildGLY(Residue originalRes) {
    }


    protected static void buildHIS(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom HIS CG C5 */
        double bond = 1.497;     /*       CH2E-C5 */
        double angle = (Math.PI / 180) * 113.8; /*  CH1E-CH2E-C5 */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);


        /* Atom HIS ND1 NH1 */
        bond = 1.378;     /*       C5-NH1 */
        angle = (Math.PI / 180) * 122.7; /*  CH2E-C5-NH1 */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "ND1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom HIS HD1 NH1 */
        bond = 1.0;     /*       C5-NH1 */
        angle = (Math.PI / 180) * 125; /*  CH2E-C5-NH1 */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "ND1"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "HD1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom HIS CE1 CRH */
        bond = 1.345;     /*     NH1-CRH */
        angle = (Math.PI / 180) * 109.0; /*  C5-NH1-CRH */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "ND1"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom HIS NE2 NR */
        bond = 1.319;     /*      CRH-NR */
        angle = (Math.PI / 180) * 111.7; /*  NH1-CRH-NR */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "CE1"), getAtom(originalRes.getAtoms(), "ND1"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "NE2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom HIS HE2 */
        bond = 1.0;     /*       C5-NH1 */
        angle = (Math.PI / 180) * 125.57; /*  CH2E-C5-NH1 */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "NE2"), getAtom(originalRes.getAtoms(), "CE1"), getAtom(originalRes.getAtoms(), "ND1"));

        getAtom(originalRes.getAtoms(), "HE2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom HIS CD2 CR1E */
        bond = 1.356;     /*      NR-CR1E */
        angle = (Math.PI / 180) * 130.69; /*  CRH-NR-CR1E */

        getAtom_xyza(xyza, bond, angle, chi[1] - Math.PI,
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildILE(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom ILE CG1 CH2E */
        double bond = 1.530;     /*       CH1E-CH2E */
        double angle = (Math.PI / 180) * 110.4; /*  CH1E-CH1E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ILE CG2 CH3E */
        bond = 1.521;     /*       CH1E-CH3E */
        angle = (Math.PI / 180) * 110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[0] - (Math.PI * 2.0 / 3.0),
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ILE CD CH3E */
        bond = 1.513;     /*       CH2E-CH3E */
        angle = (Math.PI / 180) * 113.8; /*  CH1E-CH2E-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG1"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD1").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildLYS(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom LYS CG CH2E */
        double bond = 1.520;     /*       CH2E-CH2E */
        double angle = (Math.PI / 180) * 114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom LYS CD CH2E */
        bond = 1.520;     /*       CH2E-CH2E */
        angle = (Math.PI / 180) * 111.3; /*  CH2E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom LYS CE CH2E */
        bond = 1.520;     /*       CH2E-CH2E */
        angle = (Math.PI / 180) * 111.3; /*  CH2E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[2],
                getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom LYS NZ NH3 */
        bond = 1.489;     /*       CH2E-NH3 */
        angle = (Math.PI / 180) * 111.9; /*  CH2E-CH2E-NH3 */

        getAtom_xyza(xyza, bond, angle, chi[3],
                getAtom(originalRes.getAtoms(), "CE"), getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "NZ").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildLEU(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom LEU CG CH1E */
        double bond = 1.530; /*       CH1E-CH2E */
        double angle = (Math.PI / 180) * 116.3; /*  CH1E-CH2E-CH1E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom LEU CD1 CH3E */
        bond = 1.521; /*       CH1E-CH3E */
        angle = (Math.PI / 180) * 110.7; /*  CH2E-CH1E-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom LEU CD2 CH3E */
        bond = 1.521; /*       CH1E-CH3E */
        angle = (Math.PI / 180) * 110.7; /*  CH2E-CH1E-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[1] + (Math.PI * 2.0 / 3.0),
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildMET(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom MET CG CH2E */
        double bond = 1.520;     /*       CH2E-CH2E */
        double angle = (Math.PI / 180) * 114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom MET SD SM */
        bond = 1.803;     /*       CH2E-SM */
        angle = (Math.PI / 180) * 112.7; /*  CH2E-CH2E-SM */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "SD").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom MET CE CH3E */
        bond = 1.791;     /*       SM-CH3E */
        angle = (Math.PI / 180) * 100.9; /*  CH2E-SM-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[2],
                getAtom(originalRes.getAtoms(), "SD"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildASN(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom ASN CG C */
        double bond = 1.516;     /*       CH2E-C */
        double angle = (Math.PI / 180) * 112.6; /*  CH1E-CH2E-C */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);


        /* Atom ASN OD1 O */
        bond = 1.231;     /*       C-O */
        angle = (Math.PI / 180) * 120.8; /*  CH2E-C-O */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "OD1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ASN ND2 NH2 */
        bond = 1.328;     /*       C-NH2 */
        angle = (Math.PI / 180) * 116.4; /*  CH2E-C-NH2 */

        getAtom_xyza(xyza, bond, angle, chi[1] - Math.PI,
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "ND2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ASN HD21 NH2 */
        bond = 1.0;     /*       C-NH2 */
        angle = (Math.PI / 180) * 119; /*  CH2E-C-NH2 */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "ND2"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "HD21").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ASN HD22 NH2 */
        bond = 1.0;     /*       C-NH2 */
        angle = (Math.PI / 180) * 119; /*  CH2E-C-NH2 */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "ND2"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "HD22").setXYZ(xyza[0], xyza[1], xyza[2]);


    }

    protected static void buildPRO(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom PRO CG CH2P */
        double bond = 1.492;     /*       CH2E-CH2P */
        double angle = (Math.PI / 180) * 104.5; /*  CH1E-CH2E-CH2P */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom PRO CD CH2P */
        bond = 1.503;     /*       CH2P-CH2P */
        angle = (Math.PI / 180) * 106.1; /*  CH1E-CH2P-CH2P */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildGLN(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom GLN CG CH2E */
        double bond = 1.520;     /*       CH2E-CH2E */
        double angle = (Math.PI / 180) * 114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);


        /* Atom GLU CD C */
        bond = 1.516;     /*       CH2E-C */
        angle = (Math.PI / 180) * 112.6; /*  CH2E-CH2E-C */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom GLN OE1 O */
        bond = 1.231;     /*       C-O */
        angle = (Math.PI / 180) * 120.8; /*  CH2E-C-O */

        getAtom_xyza(xyza, bond, angle, chi[2],
                getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "OE1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom GLN NE2 NH2 */
        bond = 1.328;     /*       C-NH2 */
        angle = (Math.PI / 180) * 116.4; /*  CH2E-C-NH2 */

        getAtom_xyza(xyza, bond, angle, chi[2] - Math.PI,
                getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "NE2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom GLN NE21 NH2 */
        bond = 1.0;     /*       C-NH2 */
        angle = (Math.PI / 180) * 119.0; /*  CH2E-C-NH2 */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "NE2"), getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "HE21").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom GLN NE22 NH2 */
        bond = 1.0;     /*       C-NH2 */
        angle = (Math.PI / 180) * 119.0; /*  CH2E-C-NH2 */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "NE2"), getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "HE22").setXYZ(xyza[0], xyza[1], xyza[2]);
    }


    protected static void buildARG(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom LYS CG CH2E */
        double bond = 1.520;     /*       CH2E-CH2E */
        double angle = (Math.PI / 180) * 114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);


        /* Atom LYS CD CH2E */
        bond = 1.520;     /*       CH2E-CH2E */
        angle = (Math.PI / 180) * 111.3; /*  CH2E-CH2E-CH2E */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ARG NE NH1 */
        bond = 1.461;     /*      CH2E-NH1 */
        angle = (Math.PI / 180) * 112.0; /* CH2E-CH2E-NH1 */

        getAtom_xyza(xyza, bond, angle, chi[2],
                getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "NE").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ARG HE */
        bond = 1.0;     /*      CH2E-NH1 */
        angle = (Math.PI / 180) * 117.58; /* CH2E-CH2E-NH1 */

        getAtom_xyza(xyza, bond, angle, chi[3] - Math.PI,
                getAtom(originalRes.getAtoms(), "NE"), getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "HE").setXYZ(xyza[0], xyza[1], xyza[2]);


        /* Atom ARG CZ C */
        bond = 1.329;     /*       NH1-C */
        angle = (Math.PI / 180) * 124.2; /*  CH2E-NH1-C */

        getAtom_xyza(xyza, bond, angle, chi[3],
                getAtom(originalRes.getAtoms(), "NE"), getAtom(originalRes.getAtoms(), "CD"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "CZ").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ARG NH1 NC2 */
        bond = 1.326;     /*      C-NC2 */
        angle = (Math.PI / 180) * 120.0; /*  NH1-C-NC2 */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "CZ"), getAtom(originalRes.getAtoms(), "NE"), getAtom(originalRes.getAtoms(), "CD"));

        getAtom(originalRes.getAtoms(), "NH1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom ARG NH2 NH2 */
        bond = 1.326;     /*      C-NC2 */
        angle = (Math.PI / 180) * 120.0; /*  NH1-C-NC2 */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CZ"), getAtom(originalRes.getAtoms(), "NE"), getAtom(originalRes.getAtoms(), "CD"));

        getAtom(originalRes.getAtoms(), "NH2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildSER(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom SER OG OH1 */
        double bond = 1.417;     /*       CH2E-OH1 */
        double angle = (Math.PI / 180) * 111.1; /*  CH1E-CH2E-OH1 */


        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "OG").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildTHR(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom THR OG1 OH1 */
        double bond = 1.433;     /*       CH1E-OH1 */
        double angle = (Math.PI / 180) * 109.6; /*  CH1E-CH1E-OH1 */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "OG1").setXYZ(xyza[0], xyza[1], xyza[2]);
        /* Atom THR CG2 CH3E */
        bond = 1.521;     /*       CH1E-CH3E */
        angle = (Math.PI / 180) * 110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[0] - (Math.PI * 2 / 3),
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildVAL(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom VAL CG1 CH3E */
        double bond = 1.521;     /*       CH1E-CH3E */
        double angle = (Math.PI / 180) * 110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG1").setXYZ(xyza[0], xyza[1], xyza[2]);
        /* Atom VAL CG2 CH3E */
        bond = 1.521;     /*       CH1E-CH3E */
        angle = (Math.PI / 180) * 110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza, bond, angle, chi[0] + (Math.PI * 2 / 3),
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildTRP(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom TRP CG C5W */
        double bond = 1.498;     /*       CH2E-C5W */
        double angle = (Math.PI / 180) * 113.6; /*  CH1E-CH2E-C5W */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP CD1 CR1E */
        bond = 1.365;     /*       C5W-CR1E */
        angle = (Math.PI / 180) * 126.934; /*  CH2E-C5W-CR1E */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP CD2 CW */
        bond = 1.433;     /*       C5W-CW */
        angle = (Math.PI / 180) * 126.578; /*  CH2E-C5W-CW */

        getAtom_xyza(xyza, bond, angle, chi[1] - Math.PI,
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP NE1 NH1 */
        bond = 1.374;     /*      CR1E-NH1 */
        angle = (Math.PI / 180) * 110.2; /*  C5W-CR1E-NH1 */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CD1"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "NE1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP NE1 NH1 */
        bond = 1.0;     /*      CR1E-NH1 */
        angle = (Math.PI / 180) * 125.2; /*  C5W-CR1E-NH1 */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "NE1"), getAtom(originalRes.getAtoms(), "CD1"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "HE1").setXYZ(xyza[0], xyza[1], xyza[2]);


        /* Atom TRP CE2 CW */
        bond = 1.4114;     /*      CW-CW */
        angle = (Math.PI / 180) * 107.1; /*  C5W-CW-CW */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP CE3 CR1E */
        bond = 1.401;     /*     CW-CR1E */
        angle = (Math.PI / 180) * 133.885; /*  CW-CW-CR1E */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE3").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP CZ2 CR1W */
        bond = 1.3982;     /*     CW-CR1W */
        angle = (Math.PI / 180) * 122.332; /*  CW-CW-CR1W */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CE2"), getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "CZ2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP CZ3 CR1E */
        bond = 1.391;     /*     CR1E-CR1E */
        angle = (Math.PI / 180) * 118.654; /*  CW-CR1E-CR1E */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CE3"), getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "CZ3").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TRP CH2 CR1W */
        bond = 1.368;     /*     CR1W-CR1W */
        angle = (Math.PI / 180) * 117.5; /*  CW-CR1W-CR1W */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "CZ3"), getAtom(originalRes.getAtoms(), "CE3"), getAtom(originalRes.getAtoms(), "CD2"));

        getAtom(originalRes.getAtoms(), "CH2").setXYZ(xyza[0], xyza[1], xyza[2]);
    }

    protected static void buildTYR(Residue originalRes, double[] chi) {
        buildALA(originalRes);
        /* Atom TYR CG CY */
        double bond = 1.512;     /*       CH2E-CY */
        double angle = (Math.PI / 180) * 113.9; /*  CH1E-CH2E-CY */

        getAtom_xyza(xyza, bond, angle, chi[0],
                getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"), getAtom(originalRes.getAtoms(), "N"));

        getAtom(originalRes.getAtoms(), "CG").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TYR CD1 CR1E */
        bond = 1.389;     /*       CY-CR1E */
        angle = (Math.PI / 180) * 120.8; /*  CH2E-CY-CR1E */

        getAtom_xyza(xyza, bond, angle, chi[1],
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TYR CD2 CR1E */
        bond = 1.389;     /*      CY-CR1E */
        angle = (Math.PI / 180) * 120.8; /* CH2E-CY-CR1E */

        getAtom_xyza(xyza, bond, angle, chi[1] - Math.PI,
                getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"), getAtom(originalRes.getAtoms(), "CA"));

        getAtom(originalRes.getAtoms(), "CD2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TYR CE1 CR1E */
        bond = 1.382;     /*     CR1E-CR1E */
        angle = (Math.PI / 180) * 121.2; /*  CY-CR1E-CR1E */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CD1"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE1").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TYR CE2 CR1E */
        bond = 1.382;     /*     CR1E-CR1E */
        angle = (Math.PI / 180) * 121.2; /*  CY-CR1E-CR1E */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"), getAtom(originalRes.getAtoms(), "CB"));

        getAtom(originalRes.getAtoms(), "CE2").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TYR CZ CY2 */
        bond = 1.378;     /*       CR1E-CY2 */
        angle = (Math.PI / 180) * 119.6; /*  CR1E-CR1E-CY2 */

        getAtom_xyza(xyza, bond, angle, 0,
                getAtom(originalRes.getAtoms(), "CE2"), getAtom(originalRes.getAtoms(), "CD2"), getAtom(originalRes.getAtoms(), "CG"));

        getAtom(originalRes.getAtoms(), "CZ").setXYZ(xyza[0], xyza[1], xyza[2]);

        /* Atom TYR OH OH1 */
        bond = 1.376;     /*       CY2-OH1 */
        angle = (Math.PI / 180) * 119.9; /*  CR1E-CY2-OH1 */

        getAtom_xyza(xyza, bond, angle, Math.PI,
                getAtom(originalRes.getAtoms(), "CZ"), getAtom(originalRes.getAtoms(), "CE1"), getAtom(originalRes.getAtoms(), "CD1"));

        getAtom(originalRes.getAtoms(), "OH").setXYZ(xyza[0], xyza[1], xyza[2]);
    }


    /**
     * Returns the specified atom from a residue.
     */
    private static Atom getAtom(AtomList al, String atomName) {
        Iterator atoms;
        Atom atom;
        atoms = al.iterator();
        while ((atom = (Atom) atoms.next()) != null)
            if (atom.name.equals(atomName)) {
                return atom;
            }
        return null;
    }

}

