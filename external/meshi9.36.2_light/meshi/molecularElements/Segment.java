/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.molecularElements.SegmentAttribute;
import meshi.molecularElements.Residue;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 18/01/2010
 * Time: 07:12:14
 */
public class Segment extends ArrayList<Residue> {
    private int ID;

    public int ID() {
        return ID;
    }

    public Segment(int ID) {
        this.ID = ID;
    }

    public boolean add(Residue residue) {
        residue.addAttribute(new SegmentAttribute(this));
        return super.add(residue);
    }
}
