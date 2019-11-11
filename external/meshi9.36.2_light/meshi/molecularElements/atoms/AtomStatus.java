/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.atoms;

public enum AtomStatus {
    NORMAL, FROZEN, NOWHERE, HIDDEN, IMAGE;

    public boolean active() {
        return ((this == AtomStatus.NORMAL) || (this == AtomStatus.FROZEN));
    }

    public boolean frozen() {
        return ((this == AtomStatus.FROZEN) || (this == AtomStatus.IMAGE));
    }

    public boolean nowhere() {
        return (this == AtomStatus.NOWHERE);
    }

    public boolean image() {
        return (this == AtomStatus.IMAGE);
    }

    public boolean activeOrImage() {
        return ((this == AtomStatus.IMAGE) || active());
    }
}
