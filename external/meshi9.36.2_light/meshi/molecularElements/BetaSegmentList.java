/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.parameters.SecondaryStructure;
import meshi.util.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 25/09/2010
 * Time: 22:27:28
 * To change this template use File | Settings | File Templates.
 */
public class BetaSegmentList extends SegmentList {
    private boolean printed = false;

    public BetaSegmentList(Protein protein) {
        Segment segment = null;
        boolean in = false;
        int prevNumber = -9999;
        int i = 0;
        for (Residue residue : protein.residues()) {
            if (residue.getSecondaryStructure().equals(SecondaryStructure.SHEET)) {
                if ((in == false) || (residue.number() != prevNumber + 1)) {
                    in = true;
                    segment = new Segment(i);
                    i++;
                    add(segment);
                }
                segment.add(residue);
            } else {
                in = false;
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
