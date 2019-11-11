/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.parameters.SecondaryStructure;
import meshi.util.Utils;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 18/01/2010
 * Time: 12:26:39
 * To change this template use File | Settings | File Templates.
 */
public class SegmentList extends ArrayList<Segment> {

    private boolean printed = false;

    public SegmentList() {
        super();
    }

    public SegmentList(Protein protein) {
        Segment segment = null;
        int prevNumber = -9999;
        int i = 0;
        SecondaryStructure currentSS = null;
        for (Residue residue : protein.residues()) {
            if (residue.getSecondaryStructure().equals(currentSS) && (residue.number() == prevNumber + 1)) {
                segment.add(residue);
            } else {
                segment = new Segment(i);
                segment.add(residue);
                i++;
                add(segment);
                currentSS = residue.getSecondaryStructure();
            }
            prevNumber = residue.number();
        }
        if (!printed) {
            for (Segment seg : this)
                Utils.println("segment " + seg);
            printed = true;
        }
    }
}
