/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

public enum DistanceMode {
    INFINITE(false, false, false, false, false),
    NEW(false, false, false, false, false),
    NORMAL(false, false, false, false, false),
    FREE(false, false, false, true, false),
    BONDED(false, true, false, false, false),
    FROZEN(true, false, false, false, false),
    BONDED_FROZEN(true, true, false, false, false),
    MIRROR(false, false, true, false, false),
    DEAD(false, false, true, false, true);

    public final boolean frozen;
    public final boolean bonded;
    public final boolean mirror;
    public final boolean free;
    public final boolean dead;


    private DistanceMode(boolean frozen, boolean bonded, boolean mirror, boolean free, boolean dead) {
        this.frozen = frozen;
        this.bonded = bonded;
        this.mirror = mirror;
        this.free = free;
        this.dead = dead;
    }
}


