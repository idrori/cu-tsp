/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

/*
* All the aditional atributes of Distance are type of DistanceAtribute
*/
public interface MeshiAttribute {
    final static int LENNARD_JONES_ELEMENT_ATTRIBUTE = 0;
    final static int EXCLUDED_VOLUME_ELEMENT_ATTRIBUTE = 1;
    final static int HYDROGEN_BONDS_ATTRIBUTE = 2;
    final static int CN_ATTRIBUTE = 3;
    final static int SEQUENCE_ALIGNMENT_COLUMN_ATTRIBUTE = 4;
    final static int SOLVATION_ALL_ATOM_ATTRIBUTE = 5;
    final static int SOLVATE_CA_ATTRIBUTE = 6;
    final static int RESIDUE_ATTRIBUTE = 7;
    final static int GAUSSIAN_ATTRIBUTE = 8;
    final static int GAUSSIAN_ALPHA = 9;
    final static int SECONDARY_STRUCTURE_ATTRIBUTE = 10;
    final static int DISTANCE_FROM_CATALYTIC_ATTRIBUTE = 11;
    final static int BEAUTIFY_ATTRIBUTE = 12;
    final static int LOOP_RESIDUE = 13;
    final static int ORIGINAL_ATOM = 14;
    final static int SOLVATE_ROT1_ATTRIBUTE = 15;
    final static int SIDE_CHAIN_SOLVATE_ALL_ATOM_ATTRIBUTE = 16;
    final static int SUMMA_ATTRIBUTE = 17;
    final static int RESIDUE_TORSIONS_ATTRIBUTE = 18;
    final static int ENERGY_ATTRIBUTE = 19;
    final static int SEGMENT_ATTRIBUTE = 20;
    final static int SS_PREDICTION = 21;
    final static int CONSENSUS = 22;
    final static int SS_DSSP = 23;
    final static int CONFIDENCE = 24;
    final static int CHAR = 25;
    final static int SOLVATE_ALL_ATOM_ATTRIBUTE = 26;

    int key();
}
