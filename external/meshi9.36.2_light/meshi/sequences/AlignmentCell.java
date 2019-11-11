/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.util.*;

/**
 * A container for a protein element and its number.
 */
public abstract class AlignmentCell implements Attributable {
    public final int number;
    public final Object obj;
    public final String comment;
    private AttributesRack attributes = new AttributesRack();

    public AlignmentCell(Object obj, int number) {
        this.number = number;
        this.obj = obj;
        this.comment = "UNKNOWN";
    }

    public AlignmentCell(Object obj, int number, String comment) {
        this.number = number;
        this.obj = obj;
        this.comment = comment;
    }

    public Object object() {
        return obj;
    }

    public int number() {
        return number;
    }

    public abstract boolean gap();

    public AlignmentColumn column() {
        AlignmentColumn out = new AlignmentColumn(1);
        out.add(0, this);
        return out;
    }

    //----------------------------------------------

    public void addAttribute(MeshiAttribute attribute) {
        attributes.addAttribute(attribute);
    }

    public MeshiAttribute getAttribute(int key) {
        return attributes.getAttribute(key);
    }

}
