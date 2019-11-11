/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;

import java.util.Vector;

public class LoopResultVector extends Vector<LoopResults> {

    private static final long serialVersionUID = 314159L;


    public int isLoopCloseToAnyInTheVector(LoopResults loop, double rmsCriterion) {
        for (int c = 0; c < size(); c++)
            if (loop.calcRMS(get(c)) < rmsCriterion)
                return c;
        return -1;
    }


}
