/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.util.MeshiAttribute;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 21/01/2010
 * Time: 16:01:00
 * To change this template use File | Settings | File Templates.
 */
public class SegmentAttribute implements MeshiAttribute {
    private Segment segment;

    public SegmentAttribute(Segment segment) {
        this.segment = segment;
    }

    public Segment segment() {
        return segment;
    }

    public int key() {
        return MeshiAttribute.SEGMENT_ATTRIBUTE;
    }

    public String toString() {
        return "Segment " + segment.ID();
    }
}
