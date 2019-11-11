/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

/**
 * Distances are treated differently by Solvate Energy.
 */
public enum SolvationDistanceType {
    AT_LEAST_ONE_HYDROGEN, TWO_POLARS, POLAR_HYDROPHOBIC, HYDROPHOBIC_POLAR, TWO_HYDROPHOBIC
}
