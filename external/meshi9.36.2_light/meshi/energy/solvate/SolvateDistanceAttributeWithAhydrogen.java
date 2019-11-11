/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

/**

 */
public class SolvateDistanceAttributeWithAhydrogen extends SolvateDistanceAttribute {
    public final SolvateDistanceType type() {
        return SolvateDistanceType.AT_LEAST_ONE_HYDROGEN;
    }
}
