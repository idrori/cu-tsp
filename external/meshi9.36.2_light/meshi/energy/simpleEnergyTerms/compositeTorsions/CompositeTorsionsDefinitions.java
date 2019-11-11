/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions;

import meshi.energy.simpleEnergyTerms.*;


/**
 * Definitions and constants used throughout the CompositeTorsionsEnergy
 * term.
 *
 * @author El-ad David Amir
 */
public interface CompositeTorsionsDefinitions {

    /* torsion types */
    public static final int OMG = 0;
    public static final int PHI = 1;
    public static final int PSI = 2;
    public static final int CHI_1 = 3;
    public static final int CHI_2 = 4;
    public static final int CHI_3 = 5;
    public static final int CHI_4 = 6;
    public static final int TOTAL_TORSION_ANGLES = 7;

    public static final int UNIDENTIFIED_TORSION_TYPE = -1;

    /* number of sidechain torsions for each residue type */
    public static final int NUM_SIDECHAIN_TORSIONS[] = {
            /* ALA */ 0,
            /* CYS */ 1,
            /* ASP */ 2,
            /* GLU */ 3,
            /* PHE */ 2,
            /* GLY */ 0,
            /* HIS */ 2,
            /* ILE */ 2,
            /* LYS */ 4,
            /* LEU */ 2,
            /* MET */ 3,
            /* ASN */ 2,
            /* PRO */ 2,
            /* GLN */ 3,
            /* ARG */ 4,
            /* SER */ 1,
            /* THR */ 1,
            /* VAL */ 1,
            /* TRP */ 2,
            /* TYR */ 2};

    /* the OMNI amino acid type (used in propensity) */
    public static final int OMNI = 21;

    /* secondary structure types */
    public static final int HELIX = 1;
    public static final int SHEET = 2;
    public static final int COIL = 3;
    public static final int ALL = 9;

    /* splined polynomials types */
    public static final int POLYNOMIAL_PHI_PSI = 0;
    public static final int POLYNOMIAL_PHI_PSI_CHI_1 = 1;
    public static final int POLYNOMIAL_CHI_1_CHI_2 = 2;
    public static final int POLYNOMIAL_CHI_1 = 3;
    public static final int POLYNOMIAL_CHI_1_CHI_3 = 4;
    public static final int POLYNOMIAL_CHI_1_CHI_4 = 5;

    /* splined polynomials torsion angles */
    public static final int POLYNOMIAL_PHI_PSI_TORSIONS[] = {PHI, PSI};
    public static final int POLYNOMIAL_PHI_PSI_CHI_1_TORSIONS[] = {PHI, PSI, CHI_1};
    public static final int POLYNOMIAL_CHI_1_CHI_2_TORSIONS[] = {CHI_1, CHI_2};
    public static final int POLYNOMIAL_CHI_1_TORSIONS[] = {CHI_1};
    public static final int POLYNOMIAL_CHI_1_CHI_3_TORSIONS[] = {CHI_1, CHI_3};
    public static final int POLYNOMIAL_CHI_1_CHI_4_TORSIONS[] = {CHI_1, CHI_4};

    public static final int PREPRO = 20;

}


